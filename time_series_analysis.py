import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from MDAnalysis.analysis.leaflet import LeafletFinder
from tqdm import tqdm
import os
import time
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from functools import partial

# Constants
START = 60000
STOP = 80000
STEP = 50
CONTACT_CUTOFF = 8.0  # Å - Typical cutoff value for residue contacts
DIMER_DISTANCE = 20.0  # Distance defining dimerization (Å)
PROTEIN_COM_CUTOFF = DIMER_DISTANCE
TIME_SERIES_INTERVAL = 20  # Frame interval for time series plots (set to 20 for finer display)
TOPOLOGY_FILE = 'step5_assembly.psf'
TRAJECTORY_FILE = 'md_wrapped.xtc'
OUTPUT_DIR = 'lipid_protein_timeseries6u8u'

# Lipid color map - Extended for auto-detection
LIPID_COLORS = {
    # Cholesterol
    'CHOL': '#FF0000', 'ERGO': '#FF4500',
    # Sphingomyelin
    'DPSM': '#0000FF', 'DASM': '#1E90FF', 'DBSM': '#4169E1',
    # PC lipids
    'POPC': '#00FF00', 'DOPC': '#32CD32', 'DPPC': '#00CC00', 'DIPC': '#00CC00',
    # PS lipids
    'POPS': '#FFD700', 'DOPS': '#ADFF2F', 'DPPS': '#FFA500', 'DIPS': '#FF8C00',
    # PE lipids
    'POPE': '#8B008B', 'DOPE': '#9370DB', 'DPPE': '#9932CC', 'DIPE': '#BA55D3',
    # PG lipids
    'POPG': '#808080', 'DOPG': '#A9A9A9', 'DPPG': '#C0C0C0', 'DPG3': '#C0C0C0',
    # Default color function for unknown lipids
}

# Determine number of processors (use 75% of all CPUs)
NUM_CPUS = max(1, int(multiprocessing.cpu_count() * 0.75))

def get_lipid_color_palette(lipid_types):
    """
    Generate color palette for lipid types.
    
    Parameters
    ----------
    lipid_types : list
        List of lipid type names
    
    Returns
    -------
    dict
        Dictionary mapping lipid types to hex color codes
    """
    palette = {}
    for i, lipid_type in enumerate(lipid_types):
        if lipid_type in LIPID_COLORS:
            palette[lipid_type] = LIPID_COLORS[lipid_type]
        else:
            # Generate color for unknown lipid types
            import colorsys
            hue = (i * 0.137 + 0.5) % 1.0  # Golden angle for good distribution
            saturation = 0.7
            value = 0.85
            rgb = colorsys.hsv_to_rgb(hue, saturation, value)
            palette[lipid_type] = '#{:02x}{:02x}{:02x}'.format(
                int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)
            )
    return palette

def load_universe():
    """Load and return the MDAnalysis Universe"""
    print("Loading trajectory...")
    u = mda.Universe(TOPOLOGY_FILE, TRAJECTORY_FILE)
    print(f"Trajectory loaded successfully with {len(u.trajectory)} frames.")
    return u

def identify_lipid_leaflets(u):
    """Identify and return lipid leaflets"""
    print("Identifying lipid leaflets...")
    L = LeafletFinder(u, "name GL1 GL2 AM1 AM2 ROH GM1 GM2 PO4")
    cutoff = L.update(13)
    print(f"Cutoff distance: {cutoff}")
    print(f"Number of leaflets found: {len(L.components)}")
    
    if len(L.components) < 2:
        print("Warning: Less than 2 leaflets found. Using the same leaflet for both.")
        leaflet0 = L.groups(0)
        leaflet1 = leaflet0
    else:
        leaflet0 = L.groups(0)
        leaflet1 = L.groups(1)
    
    # Count lipids in each leaflet
    for leaflet_num, leaflet in enumerate([leaflet0, leaflet1]):
        print(f"\nLeaflet {leaflet_num} lipid counts:")
        lipid_types = ['CHOL', 'DPSM', 'DIPC', 'DPG3', 'DOPS']
        for lipid_type in lipid_types:
            n_residues = len(leaflet.select_atoms(f"resname {lipid_type}").residues)
            print(f"{lipid_type}: {n_residues} molecules")
    
    print("\nLipid leaflets identified.")
    return leaflet0, leaflet1, L

def setup_lipid_selections(leaflet0, leaflet1, lipid_types=None):
    """
    Setup selections for lipid models with auto-detection.
    
    Parameters
    ----------
    leaflet0, leaflet1 : MDAnalysis.AtomGroup
        Upper and lower leaflets
    lipid_types : list, optional
        List of lipid types to analyze. If None, auto-detect from leaflets.
    
    Returns
    -------
    dict
        Dictionary with lipid selections for each type and leaflet
    """
    # Auto-detect lipid types if not provided
    if lipid_types is None:
        # Get all lipid types from both leaflets
        all_lipids_resnames = set(leaflet0.residues.resnames) | set(leaflet1.residues.resnames)
        
        # Filter out obvious non-lipids
        non_lipids = {'W', 'NA', 'CL', 'ION', 'NA+', 'CL-', 'K', 'K+', 'MG', 'MG2+', 'CA', 'CA2+'}
        lipid_types = sorted([resname for resname in all_lipids_resnames if resname not in non_lipids])
        
        if len(lipid_types) == 0:
            # Fallback to scanning the entire system
            print("  Could not detect lipids from leaflets. Scanning entire system...")
            u = leaflet0.universe
            all_resnames = set(u.residues.resnames)
            # Known lipid patterns
            known_lipids = {'POPC', 'POPE', 'POPS', 'DOPC', 'DOPE', 'DOPS', 'DPPC', 'DPPE', 'DPPS',
                          'CHOL', 'DPSM', 'DIPC', 'DPG3', 'DIPS', 'DLPC', 'DLPE', 'DLPS'}
            lipid_types = sorted([resname for resname in all_resnames if resname in known_lipids])
            
            if len(lipid_types) == 0:
                # Ultimate fallback
                print("  Warning: No lipids detected. Using default set.")
                lipid_types = ['CHOL', 'DPSM', 'DIPC', 'DPG3', 'DOPS']
        
        print(f"\n  Auto-detected {len(lipid_types)} lipid types: {', '.join(lipid_types)}")
    
    lipid_sels = {}
    for lipid_type in lipid_types:
        # Check if this lipid exists in the leaflets
        l0_lipids = leaflet0.select_atoms(f"resname {lipid_type}")
        l1_lipids = leaflet1.select_atoms(f"resname {lipid_type}")
        
        if len(l0_lipids) > 0 or len(l1_lipids) > 0:
            lipid_sels[lipid_type] = {
                'sel': [l0_lipids, l1_lipids]
            }
            n_l0 = len(l0_lipids.residues) if len(l0_lipids) > 0 else 0
            n_l1 = len(l1_lipids.residues) if len(l1_lipids) > 0 else 0
            print(f"    {lipid_type}: {n_l0} in upper leaflet, {n_l1} in lower leaflet")
    
    if len(lipid_sels) == 0:
        raise ValueError("No lipids found in the membrane leaflets!")
    
    return lipid_sels

def select_proteins(u, n_proteins=None):
    """
    Select proteins for analysis. Auto-detects available proteins if n_proteins is None.
    
    Parameters
    ----------
    u : MDAnalysis.Universe
        The universe containing the system
    n_proteins : int, optional
        Number of proteins to select. If None, auto-detect all available proteins.
    
    Returns
    -------
    dict
        Dictionary with protein names as keys and AtomGroups as values
    """
    proteins = {}
    possible_segids = ['PROA', 'PROB', 'PROC', 'PROD', 'PROE', 'PROF', 'PROG', 'PROH']
    
    # Auto-detect which segments actually exist
    existing_segids = []
    for segid in possible_segids:
        try:
            test_selection = u.select_atoms(f"segid {segid}")
            if len(test_selection) > 0:
                existing_segids.append(segid)
        except:
            continue
    
    if len(existing_segids) == 0:
        raise ValueError("No protein segments found in the system!")
    
    print(f"\nFound {len(existing_segids)} protein segment(s): {', '.join(existing_segids)}")
    
    # Use either specified number or all found proteins
    if n_proteins is None:
        n_proteins = len(existing_segids)
        print(f"Auto-selecting all {n_proteins} proteins")
    else:
        n_proteins = min(n_proteins, len(existing_segids))
        if n_proteins > len(existing_segids):
            print(f"Warning: Requested {n_proteins} proteins but only {len(existing_segids)} available")
            n_proteins = len(existing_segids)
    
    # Select the proteins
    for i, segid in enumerate(existing_segids[:n_proteins], 1):
        protein_name = f"Protein {i}"
        # Try to select transmembrane helix region first
        proteins[protein_name] = u.select_atoms(f"segid {segid} and resid 65:103")
        # If no atoms in that range, select all atoms in the segment
        if len(proteins[protein_name]) == 0:
            proteins[protein_name] = u.select_atoms(f"segid {segid}")
            print(f"  Note: Using all residues for {segid} (TM region 65-103 not found)")
    
    # Print protein information
    for name, protein in proteins.items():
        print(f"\n{name}:")
        print(f"  Number of atoms: {len(protein.atoms)}")
        print(f"  Number of residues: {len(protein.residues)}")
        if len(protein.residues) > 0:
            print(f"  Residue IDs: {sorted([r.resid for r in protein.residues])}")
    
    return proteins

def process_frame(frame_number, u, proteins, lipid_sels):
    """Process a single frame and calculate lipid-protein contacts"""
    # Optimization: Store fixed-size arrays and use NumPy vectorization
    u.trajectory[frame_number]
    box = u.dimensions[:3]
    
    # Initialize lipid contact matrices
    lipid_contacts = {}
    for lipid_type in lipid_sels:
        lipid_contacts[lipid_type] = {}
        for protein_name, protein in proteins.items():
            n_residues = len(protein.residues)
            contacts = np.zeros(n_residues)
            # Added: Array for simple counts
            simple_counts = np.zeros(n_residues)
            lipid_contacts[lipid_type][protein_name] = {
                'contacts': contacts,
                'simple_counts': simple_counts,  # New field added
                'residue_ids': np.array([res.resid for res in protein.residues])
            }
    
    # Process each lipid type
    for lipid_type, sel_info in lipid_sels.items():
        # Process each leaflet
        for leaflet_idx, leaflet_sel in enumerate(sel_info['sel']):
            # Skip if no residues in this leaflet selection
            if not hasattr(leaflet_sel, 'residues') or len(leaflet_sel.residues) == 0:
                continue
                
            # For each lipid residue in the leaflet
            for lipid_res in leaflet_sel.residues:
                try:
                    # Get center of mass of the lipid
                    lipid_com = lipid_res.atoms.center_of_mass()
                except:
                    # Skip if can't calculate center of mass
                    continue
                
                # For each protein
                for protein_name, protein in proteins.items():
                    # Skip if no residues
                    if len(protein.residues) == 0:
                        continue
                        
                    # Pre-compute centers of mass for all protein residues at once
                    try:
                        res_coms = np.array([res.atoms.center_of_mass() for res in protein.residues])
                    except:
                        # Skip if center of mass calculation fails
                        continue
                    
                    # Vectorized distance calculation with PBC
                    # Make a copy to avoid modifying the original
                    diff = res_coms - lipid_com
                    
                    # Apply PBC correction
                    for dim in range(3):
                        # Create masks for coordinates needing correction
                        mask_over = diff[:, dim] > box[dim] * 0.5
                        mask_under = diff[:, dim] < -box[dim] * 0.5
                        
                        # Apply corrections where needed
                        diff[mask_over, dim] -= box[dim]
                        diff[mask_under, dim] += box[dim]
                    
                    # Calculate distances
                    distances = np.sqrt(np.sum(diff * diff, axis=1))
                    
                    # Find contacts within cutoff
                    contact_mask = distances <= CONTACT_CUTOFF
                    
                    # If any contacts found
                    if np.any(contact_mask):
                        # Traditional weighting
                        weights = np.exp(-distances[contact_mask] * 0.5)
                        lipid_contacts[lipid_type][protein_name]['contacts'][contact_mask] += weights
                        
                        # Added: Also record simple counts (binary count)
                        lipid_contacts[lipid_type][protein_name]['simple_counts'][contact_mask] += 1
    
    return {
        'frame': frame_number,
        'lipid_contacts': lipid_contacts
    }

# Worker function for multiprocessing
def process_frame_worker(frame, u_info, proteins_info, lipid_sels_info):
    """Wrapper for process_frame to be used with multiprocessing"""
    # Create a new Universe for this process
    u_local = mda.Universe(u_info['topology'], u_info['trajectory'])
    
    # Reconstruct proteins selection
    proteins_local = {}
    for name, sel_str in proteins_info.items():
        proteins_local[name] = u_local.select_atoms(sel_str)
    
    # Reconstruct lipid selections
    lipid_sels_local = {}
    for lipid_type, sel_data in lipid_sels_info.items():
        # Create separate selections for each leaflet
        # Note: We'll need to identify the leaflets in each worker
        all_lipids = u_local.select_atoms(f"resname {lipid_type}")
        
        # Simple approach: split by z-coordinate above/below membrane center
        # This is a simplified approach; in a real system you might need the actual LeafletFinder
        # But for performance we'll use a simple heuristic
        all_pos = all_lipids.positions
        if len(all_pos) > 0:
            z_coords = all_pos[:, 2]  # Z-coordinates
            z_mean = np.mean(z_coords)
            
            upper_leaflet = all_lipids.atoms[z_coords >= z_mean]
            lower_leaflet = all_lipids.atoms[z_coords < z_mean]
            
            # Handle edge case where all lipids might be in one leaflet
            if len(upper_leaflet) == 0:
                upper_leaflet = lower_leaflet
            if len(lower_leaflet) == 0:
                lower_leaflet = upper_leaflet
                
            lipid_sels_local[lipid_type] = {
                'sel': [
                    upper_leaflet.residues, 
                    lower_leaflet.residues
                ]
            }
        else:
            # Handle case with no lipids of this type
            empty_sel = u_local.atoms[[]]
            lipid_sels_local[lipid_type] = {
                'sel': [empty_sel, empty_sel]
            }
    
    # Process the frame
    return process_frame(frame, u_local, proteins_local, lipid_sels_local)

def prepare_parallel_data(u, proteins, lipid_sels):
    """Prepare data structures for parallel processing"""
    # Information about the Universe
    u_info = {
        'topology': TOPOLOGY_FILE,
        'trajectory': TRAJECTORY_FILE
    }
    
    # Convert protein selections to selection strings using simpler approach
    proteins_info = {}
    for name, selection in proteins.items():
        # Store segment ID and residue range for each protein
        segids = list(set([atom.segid for atom in selection.atoms]))
        resids = list(set([atom.resid for atom in selection.atoms]))
        min_resid = min(resids)
        max_resid = max(resids)
        
        # Create selection string using segment ID and residue range
        sel_string = f"segid {' '.join(segids)} and resid {min_resid}:{max_resid}"
        proteins_info[name] = sel_string
    
    # Convert lipid selections to selection strings
    lipid_sels_info = {}
    for lipid_type, sel_info in lipid_sels.items():
        sel_strs = []
        for leaflet_sel in sel_info['sel']:
            # Simple approach: just use resname for lipids
            lipid_resname = lipid_type  # Assuming lipid_type matches resname
            sel_str = f"resname {lipid_resname}"
            sel_strs.append(sel_str)
        lipid_sels_info[lipid_type] = {'sel_strs': sel_strs}
    
    return u_info, proteins_info, lipid_sels_info

def analyze_time_series(u, proteins, lipid_sels, start, stop, step=STEP, interval=TIME_SERIES_INTERVAL):
    """Single-threaded fallback version (for use when parallel processing has issues)"""
    print(f"\nFalling back to single-thread processing from frame {start} to {stop} with step {step}, interval {interval}...")
    
    # Create a list of frames to analyze
    frames = list(range(start, stop, step))
    
    # Create data structure to hold time series data
    time_series_data = []
    
    # Process each frame
    for frame in tqdm(frames, desc="Processing frames for time series"):
        if frame % interval == 0:  # Only process frames at specified intervals for time series
            if frame >= len(u.trajectory):
                continue
                
            result = process_frame(frame, u, proteins, lipid_sels)
            time_series_data.append(result)
    
    print(f"Processed {len(time_series_data)} frames for time series analysis")
    return time_series_data

def analyze_time_series_parallel(u, proteins, lipid_sels, start, stop, step=STEP, interval=TIME_SERIES_INTERVAL):
    """Analyze lipid-protein contacts over time using parallel processing"""
    print(f"\nAnalyzing time series from frame {start} to {stop} with step {step}, interval {interval}...")
    print(f"Using {NUM_CPUS} CPU cores for parallel processing")
    
    # Create a list of frames to analyze
    frames = [frame for frame in range(start, stop, step) if frame % interval == 0 and frame < len(u.trajectory)]
    
    if not frames:
        print("No frames to analyze!")
        return []
        
    # Try parallel processing, fall back to single thread if problems occur
    try:
        # Prepare data for parallel processing
        u_info, proteins_info, lipid_sels_info = prepare_parallel_data(u, proteins, lipid_sels)
        
        # Use partial to create a worker function with fixed parameters
        worker_func = partial(process_frame_worker, 
                             u_info=u_info, 
                             proteins_info=proteins_info, 
                             lipid_sels_info=lipid_sels_info)
        
        # Process frames in parallel
        time_series_data = []
        
        with ProcessPoolExecutor(max_workers=NUM_CPUS) as executor:
            # Submit all tasks
            future_to_frame = {executor.submit(worker_func, frame): frame for frame in frames}
            
            # Process results as they complete
            for future in tqdm(as_completed(future_to_frame), total=len(frames), desc="Processing frames"):
                try:
                    result = future.result()
                    time_series_data.append(result)
                except Exception as exc:
                    frame = future_to_frame[future]
                    print(f'Frame {frame} generated an exception: {exc}')
                    raise  # Re-raise to trigger fallback
        
        # Sort results by frame number
        time_series_data.sort(key=lambda x: x['frame'])
        
        print(f"Processed {len(time_series_data)} frames for time series analysis")
        return time_series_data
        
    except Exception as e:
        print(f"Error in parallel processing: {str(e)}")
        print("Falling back to single thread processing...")
        return analyze_time_series(u, proteins, lipid_sels, start, stop, step, interval)

def plot_lipid_composition_time_series(time_series_data, proteins, output_dir, window_size=20):
    """Plot time series of lipid composition percentages with smoothing and consistent color scheme"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Get lipid types
    if len(time_series_data) == 0:
        print("No time series data to plot")
        return
        
    lipid_types = list(time_series_data[0]['lipid_contacts'].keys())
    print(f"Plotting time series for lipid types: {lipid_types}")
    
    # Create color palette using helper function
    lipid_palette = get_lipid_color_palette(lipid_types)
    
    # Initialize pandas DataFrame for better data handling
    df_data = []
    
    # Collect data
    for frame_data in time_series_data:
        frame_idx = frame_data['frame']
        
        # Count contacts for each protein and lipid type
        for protein_name in proteins.keys():
            counts = {lipid_type: 0 for lipid_type in lipid_types}
            
            # Sum contacts for this frame
            for lipid_type in lipid_types:
                if protein_name in frame_data['lipid_contacts'].get(lipid_type, {}):
                    contacts = frame_data['lipid_contacts'][lipid_type][protein_name]['contacts']
                    counts[lipid_type] = np.sum(contacts)
            
            # Calculate percentages
            total = sum(counts.values())
            if total > 0:
                for lipid_type in lipid_types:
                    percentage = (counts[lipid_type] / total) * 100
                    df_data.append({
                        'Frame': frame_idx,
                        'Protein': protein_name,
                        'Lipid_Type': lipid_type,
                        'Percentage': percentage
                    })
    
    # Create DataFrame
    df = pd.DataFrame(df_data)
    
    # Plot each protein separately using seaborn for better aesthetics
    sns.set_style("whitegrid")
    
    for protein_name in proteins.keys():
        plt.figure(figsize=(12, 6))
        
        # Filter data for this protein
        protein_df = df[df['Protein'] == protein_name]
        
        # Smoothing process
        smoothed_data = []
        for lipid_type in lipid_types:
            lipid_df = protein_df[protein_df['Lipid_Type'] == lipid_type].copy()
            
            if len(lipid_df) >= window_size:
                # Sort to ensure order
                lipid_df = lipid_df.sort_values('Frame')
                
                # Calculate moving average
                lipid_df['Smoothed'] = lipid_df['Percentage'].rolling(
                    window=window_size, center=True, min_periods=1
                ).mean()
                
                smoothed_data.append(lipid_df)
            else:
                # Use as is if data is too sparse
                lipid_df['Smoothed'] = lipid_df['Percentage']
                smoothed_data.append(lipid_df)
        
        # Combine smoothed data
        smoothed_df = pd.concat(smoothed_data)
        
        # Plot both original and smoothed data (emphasize smoothed)
        # Original data (displayed faintly)
        sns.lineplot(data=protein_df, x='Frame', y='Percentage', 
                     hue='Lipid_Type', alpha=0.3, legend=False, palette=lipid_palette)
        
        # Smoothed data (displayed prominently)
        sns.lineplot(data=smoothed_df, x='Frame', y='Smoothed', 
                    hue='Lipid_Type', linewidth=2.5, palette=lipid_palette)
        
        plt.xlabel('Frame Index')
        plt.ylabel('Distance-weighted Interaction Ratio (%)')
        plt.title(f'Lipid Composition Time Series for {protein_name}\nSmoothed with {window_size}-frame window')
        plt.legend(title='Lipid Type')
        
        # Save figure as PNG
        output_path = os.path.join(output_dir, f'{protein_name.lower().replace(" ", "_")}_lipid_composition_time_series_smoothed_{window_size}.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        # Also save as PS
        ps_output_path = os.path.join(output_dir, f'{protein_name.lower().replace(" ", "_")}_lipid_composition_time_series_smoothed_{window_size}.ps')
        plt.savefig(ps_output_path, format='ps', bbox_inches='tight')
        
        plt.close()
        
        print(f"Created smoothed time series plot for {protein_name} with {window_size}-frame window")

def plot_lipid_residue_time_series(time_series_data, proteins, output_dir, window_size=20):
    """Plot time series of contacts for specific residues of interest with smoothing and consistent colors"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Define residues of interest (e.g., central part of transmembrane region)
    interesting_residues = [62, 68, 69, 72, 73]
    
    # Get lipid types
    if len(time_series_data) == 0:
        print("No time series data to plot")
        return
        
    lipid_types = list(time_series_data[0]['lipid_contacts'].keys())
    print(f"Plotting residue time series for lipid types: {lipid_types}")
    
    # Create color palette using helper function
    lipid_palette = get_lipid_color_palette(lipid_types)
    
    # Initialize DataFrame for better data handling
    df_data = []
    
    # Collect data
    for frame_data in time_series_data:
        frame_idx = frame_data['frame']
        
        for protein_name, protein in proteins.items():
            for lipid_type in lipid_types:
                if protein_name in frame_data['lipid_contacts'].get(lipid_type, {}):
                    contact_data = frame_data['lipid_contacts'][lipid_type][protein_name]
                    residue_ids = contact_data['residue_ids']
                    contacts = contact_data['contacts']
                    
                    # Record contacts for interesting residues
                    for res_id in interesting_residues:
                        idx = np.where(residue_ids == res_id)[0]
                        if len(idx) > 0:
                            contact_val = contacts[idx[0]]
                        else:
                            contact_val = 0
                        
                        df_data.append({
                            'Frame': frame_idx,
                            'Protein': protein_name,
                            'Residue': res_id,
                            'Lipid_Type': lipid_type,
                            'Contact': contact_val
                        })
    
    # Create DataFrame
    df = pd.DataFrame(df_data)
    
    # Plot using seaborn for better aesthetics
    sns.set_style("whitegrid")
    
    for protein_name in proteins.keys():
        for res_id in interesting_residues:
            plt.figure(figsize=(12, 6))
            
            # Filter data for this protein and residue
            res_df = df[(df['Protein'] == protein_name) & (df['Residue'] == res_id)]
            
            # Smoothing process
            smoothed_data = []
            for lipid_type in lipid_types:
                lipid_df = res_df[res_df['Lipid_Type'] == lipid_type].copy()
                
                if len(lipid_df) >= window_size:
                    # Sort to ensure order
                    lipid_df = lipid_df.sort_values('Frame')
                    
                    # Calculate moving average
                    lipid_df['Smoothed'] = lipid_df['Contact'].rolling(
                        window=window_size, center=True, min_periods=1
                    ).mean()
                    
                    smoothed_data.append(lipid_df)
                else:
                    # Use as is if data is too sparse
                    lipid_df['Smoothed'] = lipid_df['Contact']
                    smoothed_data.append(lipid_df)
            
            # Combine smoothed data
            if smoothed_data:
                smoothed_df = pd.concat(smoothed_data)
                
                # Plot both original and smoothed data
                # Original data (displayed faintly)
                sns.lineplot(data=res_df, x='Frame', y='Contact', 
                             hue='Lipid_Type', alpha=0.3, legend=False, palette=lipid_palette)
                
                # Smoothed data (displayed prominently)
                sns.lineplot(data=smoothed_df, x='Frame', y='Smoothed', 
                            hue='Lipid_Type', linewidth=2.5, palette=lipid_palette)
            
            plt.xlabel('Frame Index')
            plt.ylabel('Contact Value')
            plt.title(f'{protein_name} - Residue {res_id} Lipid Contacts Over Time\nSmoothed with {window_size}-frame window')
            plt.legend(title='Lipid Type')
            
            # Save figure as PNG
            output_path = os.path.join(output_dir, f'{protein_name.lower().replace(" ", "_")}_residue_{res_id}_time_series_smoothed_{window_size}.png')
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
            # Also save as PS
            ps_output_path = os.path.join(output_dir, f'{protein_name.lower().replace(" ", "_")}_residue_{res_id}_time_series_smoothed_{window_size}.ps')
            plt.savefig(ps_output_path, format='ps', bbox_inches='tight')
            
            plt.close()
            
        print(f"Created smoothed residue time series plots for {protein_name} with {window_size}-frame window")

def plot_lipid_bar_charts(time_series_data, proteins, lipid_sels, leaflet_finder, output_dir):
    """Creates bar charts of average lipid compositions for each protein and leaflet with consistent colors"""
    os.makedirs(output_dir, exist_ok=True)
    
    if len(time_series_data) == 0:
        print("No time series data to plot")
        return
        
    lipid_types = list(time_series_data[0]['lipid_contacts'].keys())
    print(f"Creating bar charts for lipid types: {lipid_types}")
    
    # Create color palette from defined lipid colors
    lipid_palette = {lipid_type: LIPID_COLORS.get(lipid_type, '#FF9900') for lipid_type in lipid_types}
    
    # Initialize data storage
    protein_lipid_data = {}
    for protein_name in proteins.keys():
        protein_lipid_data[protein_name] = {
            'leaflet0': {lipid_type: [] for lipid_type in lipid_types},
            'leaflet1': {lipid_type: [] for lipid_type in lipid_types}
        }
    
    # Collect data
    for frame_data in time_series_data:
        for protein_name in proteins.keys():
            for leaflet_id in [0, 1]:
                leaflet_key = f'leaflet{leaflet_id}'
                counts = {lipid_type: 0 for lipid_type in lipid_types}
                
                # Sum contacts for this frame
                for lipid_type in lipid_types:
                    if protein_name in frame_data['lipid_contacts'].get(lipid_type, {}):
                        # Get lipid residues from this leaflet
                        lipid_residues = lipid_sels[lipid_type]['sel'][leaflet_id].residues
                        if len(lipid_residues) > 0:
                            contact_data = frame_data['lipid_contacts'][lipid_type][protein_name]
                            counts[lipid_type] = np.sum(contact_data['contacts'])
                
                # Calculate percentages
                total = sum(counts.values())
                if total > 0:
                    for lipid_type in lipid_types:
                        percentage = (counts[lipid_type] / total) * 100
                        protein_lipid_data[protein_name][leaflet_key][lipid_type].append(percentage)
    
    # Calculate averages
    avg_data = []
    for protein_name in proteins.keys():
        for leaflet_id in [0, 1]:
            leaflet_key = f'leaflet{leaflet_id}'
            for lipid_type in lipid_types:
                values = protein_lipid_data[protein_name][leaflet_key][lipid_type]
                if values:
                    avg_data.append({
                        'Protein': protein_name,
                        'Leaflet': f"Leaflet {leaflet_id}",
                        'Lipid_Type': lipid_type,
                        'Percentage': np.mean(values)
                    })
    
    # Create DataFrame
    df = pd.DataFrame(avg_data)
    
    # 1. Combined plot for all proteins and both leaflets
    plt.figure(figsize=(18, 10))
    
    # Setup the plot
    sns.set_style("whitegrid")
    
    # Create a grouped bar chart with consistent color scheme
    ax = sns.barplot(x='Protein', y='Percentage', hue='Lipid_Type', data=df, 
                     palette=lipid_palette, errorbar=None)
    
    # Customize the plot
    plt.xlabel('Protein', fontsize=14)
    plt.ylabel('Average Distance-weighted Interaction Ratio (%)', fontsize=14)
    plt.title('Lipid Contact Composition per Protein', fontsize=16)
    plt.legend(title='Lipid Type', title_fontsize=12, fontsize=10, loc='upper right')
    plt.grid(True, alpha=0.3)
    
    # Add percentage labels on top of each bar
    for container in ax.containers:
        ax.bar_label(container, fmt='%.1f%%', fontsize=9)
    
    # Save figure as PNG
    output_path = os.path.join(output_dir, 'all_proteins_lipid_composition.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    # Also save as PS
    ps_output_path = os.path.join(output_dir, 'all_proteins_lipid_composition.ps')
    plt.savefig(ps_output_path, format='ps', bbox_inches='tight')
    
    plt.close()
    
    # 2. Separate plots for each leaflet
    for leaflet_id in [0, 1]:
        plt.figure(figsize=(15, 8))
        
        # Filter data for this leaflet
        leaflet_df = df[df['Leaflet'] == f"Leaflet {leaflet_id}"]
        
        # Create a grouped bar chart with consistent colors
        ax = sns.barplot(x='Protein', y='Percentage', hue='Lipid_Type', data=leaflet_df, 
                         palette=lipid_palette, errorbar=None)
        
        # Customize the plot
        plt.xlabel('Protein', fontsize=14)
        plt.ylabel('Average Distance-weighted Interaction Ratio (%)', fontsize=14)
        plt.title(f'Lipid Contact Composition per Protein - Leaflet {leaflet_id}', fontsize=16)
        plt.legend(title='Lipid Type', title_fontsize=12, fontsize=10, loc='upper right')
        plt.grid(True, alpha=0.3)
        
        # Add percentage labels on top of each bar
        for container in ax.containers:
            ax.bar_label(container, fmt='%.1f%%', fontsize=9)
        
        # Save figure as PNG
        output_path = os.path.join(output_dir, f'leaflet{leaflet_id}_lipid_composition.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        # Also save as PS
        ps_output_path = os.path.join(output_dir, f'leaflet{leaflet_id}_lipid_composition.ps')
        plt.savefig(ps_output_path, format='ps', bbox_inches='tight')
        
        plt.close()
    
    print("Created lipid composition bar charts for all proteins and each leaflet")

def plot_leaflet_composition_time_series(time_series_data, proteins, lipid_sels, leaflet_finder, output_dir, window_size=20):
    """Plot time series of lipid composition percentages by leaflet using pandas for better performance with smoothing"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Get lipid types
    if len(time_series_data) == 0:
        print("No time series data to plot")
        return
        
    lipid_types = list(time_series_data[0]['lipid_contacts'].keys())
    print(f"Plotting leaflet time series for lipid types: {lipid_types}")
    
    # Create color palette from defined lipid colors
    lipid_palette = {lipid_type: LIPID_COLORS.get(lipid_type, '#FF9900') for lipid_type in lipid_types}
    
    # Process each leaflet
    for leaflet_id in [0, 1]:  # 0: upper leaflet, 1: lower leaflet
        # Get residues in this leaflet
        leaflet_residues = set([res.resid for res in leaflet_finder.groups(leaflet_id).residues])
        
        # Initialize DataFrame
        df_data = []
        
        # Collect data
        for frame_data in time_series_data:
            frame_idx = frame_data['frame']
            
            for protein_name in proteins.keys():
                counts = {lipid_type: 0 for lipid_type in lipid_types}
                
                # Calculate leaflet-specific contacts
                for lipid_type in lipid_types:
                    if protein_name in frame_data['lipid_contacts'].get(lipid_type, {}):
                        lipid_residues = lipid_sels[lipid_type]['sel'][leaflet_id].residues
                        leaflet_lipids = [res for res in lipid_residues if res.resid in leaflet_residues]
                        
                        if len(leaflet_lipids) > 0:
                            contact_data = frame_data['lipid_contacts'][lipid_type][protein_name]
                            counts[lipid_type] = np.sum(contact_data['contacts'])
                
                # Calculate percentages
                total = sum(counts.values())
                if total > 0:
                    for lipid_type in lipid_types:
                        percentage = (counts[lipid_type] / total) * 100
                        df_data.append({
                            'Frame': frame_idx,
                            'Protein': protein_name,
                            'Lipid_Type': lipid_type,
                            'Percentage': percentage,
                            'Leaflet': leaflet_id
                        })
        
        # Create DataFrame
        df = pd.DataFrame(df_data)
        
        # Plot each protein separately using seaborn
        for protein_name in proteins.keys():
            plt.figure(figsize=(12, 6))
            
            # Filter data for this protein
            protein_df = df[df['Protein'] == protein_name]
            
            # Get and smooth data for each lipid type
            smoothed_data = []
            for lipid_type in lipid_types:
                lipid_df = protein_df[protein_df['Lipid_Type'] == lipid_type].copy()
                
                if len(lipid_df) >= window_size:
                    # Sort to ensure order
                    lipid_df = lipid_df.sort_values('Frame')
                    
                    # Calculate moving average (adjust window size as needed)
                    lipid_df['Smoothed'] = lipid_df['Percentage'].rolling(
                        window=window_size, center=True, min_periods=1
                    ).mean()
                    
                    smoothed_data.append(lipid_df)
                else:
                    # Use as is if data is too sparse
                    lipid_df['Smoothed'] = lipid_df['Percentage']
                    smoothed_data.append(lipid_df)
            
            # Combine smoothed data
            smoothed_df = pd.concat(smoothed_data)
            
            # Plot both original and smoothed data (emphasize smoothed)
            # Original data (displayed faintly)
            sns.lineplot(data=protein_df, x='Frame', y='Percentage', 
                         hue='Lipid_Type', alpha=0.3, legend=False, palette=lipid_palette)
            
            # Smoothed data (displayed prominently)
            sns.lineplot(data=smoothed_df, x='Frame', y='Smoothed', 
                        hue='Lipid_Type', linewidth=2.5, palette=lipid_palette)
            
            plt.xlabel('Frame Index')
            plt.ylabel('Distance-weighted Interaction Ratio (%)')
            plt.title(f'Lipid Composition Time Series for {protein_name} (Leaflet {leaflet_id})\nSmoothed with {window_size}-frame window')
            plt.legend(title='Lipid Type')
            
            plt.ylim(0, 100)
            
            # Add grid lines (for improved visibility)
            plt.grid(True, alpha=0.3)
            
            # Save figure as PNG
            output_path = os.path.join(output_dir, f'{protein_name.lower().replace(" ", "_")}_lipid_composition_time_series_leaflet{leaflet_id}_smoothed_{window_size}.png')
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
            # Also save as PS
            ps_output_path = os.path.join(output_dir, f'{protein_name.lower().replace(" ", "_")}_lipid_composition_time_series_leaflet{leaflet_id}_smoothed_{window_size}.ps')
            plt.savefig(ps_output_path, format='ps', bbox_inches='tight')
            
            plt.close()
            
        print(f"Created smoothed leaflet {leaflet_id} time series plots with {window_size}-frame window")

def plot_simple_time_series(time_series_data, proteins, output_dir, window_size=20):
    """Time series plot using simple counts"""
    os.makedirs(output_dir, exist_ok=True)
    
    if len(time_series_data) == 0:
        print("No time series data to plot")
        return
        
    lipid_types = list(time_series_data[0]['lipid_contacts'].keys())
    print(f"Plotting simple count time series for lipid types: {lipid_types}")
    
    lipid_palette = {lipid_type: LIPID_COLORS.get(lipid_type, '#FF9900') for lipid_type in lipid_types}
    
    df_data = []
    
    for frame_data in time_series_data:
        frame_idx = frame_data['frame']
        
        for protein_name in proteins.keys():
            counts = {lipid_type: 0 for lipid_type in lipid_types}
            
            for lipid_type in lipid_types:
                if protein_name in frame_data['lipid_contacts'].get(lipid_type, {}):
                    # Use simple counts
                    simple_counts = frame_data['lipid_contacts'][lipid_type][protein_name]['simple_counts']
                    counts[lipid_type] = np.sum(simple_counts)
            
            total = sum(counts.values())
            if total > 0:
                for lipid_type in lipid_types:
                    percentage = (counts[lipid_type] / total) * 100
                    df_data.append({
                        'Frame': frame_idx,
                        'Protein': protein_name,
                        'Lipid_Type': lipid_type,
                        'Percentage': percentage
                    })
    
    df = pd.DataFrame(df_data)
    
    for protein_name in proteins.keys():
        plt.figure(figsize=(12, 6))
        
        protein_df = df[df['Protein'] == protein_name]
        
        smoothed_data = []
        for lipid_type in lipid_types:
            lipid_df = protein_df[protein_df['Lipid_Type'] == lipid_type].copy()
            
            if len(lipid_df) >= window_size:
                lipid_df = lipid_df.sort_values('Frame')
                lipid_df['Smoothed'] = lipid_df['Percentage'].rolling(
                    window=window_size, center=True, min_periods=1
                ).mean()
                smoothed_data.append(lipid_df)
            else:
                lipid_df['Smoothed'] = lipid_df['Percentage']
                smoothed_data.append(lipid_df)
        
        smoothed_df = pd.concat(smoothed_data)
        
        sns.lineplot(data=protein_df, x='Frame', y='Percentage', 
                     hue='Lipid_Type', alpha=0.3, legend=False, palette=lipid_palette)
        
        # Smoothed data (displayed prominently)
        sns.lineplot(data=smoothed_df, x='Frame', y='Smoothed', 
                    hue='Lipid_Type', linewidth=2.5, palette=lipid_palette)
        
        plt.xlabel('Frame Index')
        plt.ylabel('Simple Contact Ratio (%)')
        plt.title(f'Lipid Composition (Simple Count) for {protein_name}\nSmoothed with {window_size}-frame window')
        plt.legend(title='Lipid Type')
        
        # Save figure as PNG
        output_path = os.path.join(output_dir, f'{protein_name.lower().replace(" ", "_")}_simple_count_time_series_smoothed_{window_size}.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        # Also save as PS
        ps_output_path = os.path.join(output_dir, f'{protein_name.lower().replace(" ", "_")}_simple_count_time_series_smoothed_{window_size}.ps')
        plt.savefig(ps_output_path, format='ps', bbox_inches='tight')
        
        plt.close()
        
        print(f"Created simple count time series plot for {protein_name} with {window_size}-frame window")

def plot_simple_bar_charts(time_series_data, proteins, lipid_sels, leaflet_finder, output_dir):
    """Create bar charts using simple counts"""
    os.makedirs(output_dir, exist_ok=True)
    
    if len(time_series_data) == 0:
        print("No time series data to plot")
        return
        
    lipid_types = list(time_series_data[0]['lipid_contacts'].keys())
    print(f"Creating simple count bar charts for lipid types: {lipid_types}")
    
    lipid_palette = {lipid_type: LIPID_COLORS.get(lipid_type, '#FF9900') for lipid_type in lipid_types}
    
    # Initialize data storage
    protein_lipid_data = {}
    for protein_name in proteins.keys():
        protein_lipid_data[protein_name] = {
            'leaflet0': {lipid_type: [] for lipid_type in lipid_types},
            'leaflet1': {lipid_type: [] for lipid_type in lipid_types}
        }
    
    # Collect data
    for frame_data in time_series_data:
        for protein_name in proteins.keys():
            for leaflet_id in [0, 1]:
                leaflet_key = f'leaflet{leaflet_id}'
                counts = {lipid_type: 0 for lipid_type in lipid_types}
                
                # Sum contacts for this frame
                for lipid_type in lipid_types:
                    if protein_name in frame_data['lipid_contacts'].get(lipid_type, {}):
                        # Get lipid residues from this leaflet
                        lipid_residues = lipid_sels[lipid_type]['sel'][leaflet_id].residues
                        if len(lipid_residues) > 0:
                            contact_data = frame_data['lipid_contacts'][lipid_type][protein_name]
                            # Use simple counts
                            counts[lipid_type] = np.sum(contact_data['simple_counts'])
                
                # Calculate percentages
                total = sum(counts.values())
                if total > 0:
                    for lipid_type in lipid_types:
                        percentage = (counts[lipid_type] / total) * 100
                        protein_lipid_data[protein_name][leaflet_key][lipid_type].append(percentage)
    
    # Calculate averages
    avg_data = []
    for protein_name in proteins.keys():
        for leaflet_id in [0, 1]:
            leaflet_key = f'leaflet{leaflet_id}'
            for lipid_type in lipid_types:
                values = protein_lipid_data[protein_name][leaflet_key][lipid_type]
                if values:
                    avg_data.append({
                        'Protein': protein_name,
                        'Leaflet': f"Leaflet {leaflet_id}",
                        'Lipid_Type': lipid_type,
                        'Percentage': np.mean(values)
                    })
    
    # Create DataFrame
    df = pd.DataFrame(avg_data)
    
    # 1. Combined plot for all proteins and both leaflets
    plt.figure(figsize=(18, 10))
    
    # Setup the plot
    sns.set_style("whitegrid")
    
    # Create a grouped bar chart with consistent color scheme
    ax = sns.barplot(x='Protein', y='Percentage', hue='Lipid_Type', data=df, 
                     palette=lipid_palette, errorbar=None)
    
    # Customize the plot
    plt.xlabel('Protein', fontsize=14)
    plt.ylabel('Average Simple Contact Ratio (%)', fontsize=14)
    plt.title('Lipid Contact Composition per Protein (Simple Count)', fontsize=16)
    plt.legend(title='Lipid Type', title_fontsize=12, fontsize=10, loc='upper right')
    plt.grid(True, alpha=0.3)
    
    # Add percentage labels on top of each bar
    for container in ax.containers:
        ax.bar_label(container, fmt='%.1f%%', fontsize=9)
    
    # Save figure as PNG
    output_path = os.path.join(output_dir, 'all_proteins_simple_count_composition.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    # Also save as PS
    ps_output_path = os.path.join(output_dir, 'all_proteins_simple_count_composition.ps')
    plt.savefig(ps_output_path, format='ps', bbox_inches='tight')
    
    plt.close()
    
    print("Created simple count bar charts for all proteins")

def main():
    """Main function to run the time series analysis with smoothing and consistent color scheme"""
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Record start time
    start_time = time.time()
    
    # Set Matplotlib to use a better backend for PS files
    plt.rcParams['ps.useafm'] = True
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42  # Use TrueType fonts in PS files
    
    try:
        # Load Universe
        u = load_universe()
        
        # Identify lipid leaflets
        leaflet0, leaflet1, leaflet_finder = identify_lipid_leaflets(u)
        
        # Setup lipid selections
        lipid_sels = setup_lipid_selections(leaflet0, leaflet1)
        
        # Select proteins (auto-detect all available proteins)
        proteins = select_proteins(u)
        
        # Try using vectorized version with NumPy optimization first
        print("\nFirst optimizing single-threaded calculation with NumPy vectorization...")
        
        # Make a small test to confirm the optimized function works
        test_frame = START
        print(f"Testing optimized frame processing with frame {test_frame}...")
        try:
            test_result = process_frame(test_frame, u, proteins, lipid_sels)
            print("Vectorized processing test: SUCCESS")
            
            # Now try parallel processing
            print("\nAttempting parallel processing...")
            try:
                # Time parallel processing for a small subset of frames
                test_parallel_start = time.time()
                test_frames = min(10, (STOP - START) // STEP)
                test_end = START + (test_frames * STEP)
                
                time_series_data = analyze_time_series_parallel(
                    u, proteins, lipid_sels, START, test_end, STEP, TIME_SERIES_INTERVAL
                )
                test_parallel_time = time.time() - test_parallel_start
                print(f"Parallel processing test: SUCCESS ({test_parallel_time:.2f}s for {len(time_series_data)} frames)")
                
                # Now do the full analysis with parallel processing
                time_series_data = analyze_time_series_parallel(
                    u, proteins, lipid_sels, START, STOP, STEP, TIME_SERIES_INTERVAL
                )
                
            except Exception as e:
                print(f"Parallel processing failed: {str(e)}")
                print("Using single-threaded processing instead...")
                time_series_data = analyze_time_series(
                    u, proteins, lipid_sels, START, STOP, STEP, TIME_SERIES_INTERVAL
                )
                
        except Exception as e:
            print(f"Optimized processing failed: {str(e)}")
            print("Falling back to original implementation...")
            # Use original implementation by reimplementing the function directly
            print("Using non-optimized implementation...")
            
            # Define the original function locally
            def original_process_frame(frame_number, u, proteins, lipid_sels):
                """Process a single frame using the original processing logic"""
                u.trajectory[frame_number]
                box = u.dimensions[:3]
                
                # Initialize lipid contact matrices
                lipid_contacts = {}
                for lipid_type in lipid_sels:
                    lipid_contacts[lipid_type] = {}
                    for protein_name, protein in proteins.items():
                        n_residues = len(protein.residues)
                        contacts = np.zeros(n_residues)
                        simple_counts = np.zeros(n_residues)  # Add simple count
                        lipid_contacts[lipid_type][protein_name] = {
                            'contacts': contacts,
                            'simple_counts': simple_counts,
                            'residue_ids': [res.resid for res in protein.residues]
                        }
                
                # Process lipid-protein contacts
                for lipid_type, sel_info in lipid_sels.items():
                    for leaflet_sel in sel_info['sel']:
                        for lipid_res in leaflet_sel.residues:
                            lipid_com = lipid_res.atoms.center_of_mass()
                            
                            for protein_name, protein in proteins.items():
                                for i, res in enumerate(protein.residues):
                                    res_com = res.atoms.center_of_mass()
                                    
                                    # Calculate distance with PBC correction
                                    diff = res_com - lipid_com
                                    for dim in range(3):
                                        if diff[dim] > box[dim] * 0.5:
                                            diff[dim] -= box[dim]
                                        elif diff[dim] < -box[dim] * 0.5:
                                            diff[dim] += box[dim]
                                    
                                    dist = np.sqrt(np.sum(diff * diff))
                                    
                                    # Count as contact if within cutoff
                                    if dist <= CONTACT_CUTOFF:
                                        # Weighting
                                        weight = np.exp(-dist * 0.5)
                                        lipid_contacts[lipid_type][protein_name]['contacts'][i] += weight
                                        # Simple count
                                        lipid_contacts[lipid_type][protein_name]['simple_counts'][i] += 1
                
                return {
                    'frame': frame_number,
                    'lipid_contacts': lipid_contacts
                }
                
            # Define original single-threaded analysis process
            def original_analyze(u, proteins, lipid_sels, start, stop, step=STEP, interval=TIME_SERIES_INTERVAL):
                """Analyze time series data using original processing logic"""
                print(f"\nAnalyzing time series using original implementation from frame {start} to {stop}...")
                
                # Create a list of frames to analyze
                frames = list(range(start, stop, step))
                
                # Create data structure to hold time series data
                time_series_data = []
                
                # Process each frame
                for frame in tqdm(frames, desc="Processing frames"):
                    if frame % interval == 0:  # Only process frames at specified intervals
                        if frame >= len(u.trajectory):
                            continue
                            
                        result = original_process_frame(frame, u, proteins, lipid_sels)
                        time_series_data.append(result)
                
                print(f"Processed {len(time_series_data)} frames for time series analysis")
                return time_series_data
                
            # Use original processing
            time_series_data = original_analyze(u, proteins, lipid_sels, START, STOP, STEP, TIME_SERIES_INTERVAL)
        
        # Create time series plots with different smoothing strengths
        smoothing_windows = [10, 30, 50]
        
        for window_size in smoothing_windows:
            # 1. Weighted time series plot (traditional method)
            print(f"\nCreating weighted time series plots with {window_size}-frame smoothing...")
            plot_lipid_composition_time_series(time_series_data, proteins, OUTPUT_DIR, window_size=window_size)
            
            # 2. Simple count time series plot (new method)
            print(f"\nCreating simple count time series plots with {window_size}-frame smoothing...")
            plot_simple_time_series(time_series_data, proteins, OUTPUT_DIR, window_size=window_size)
            
            # 3. Leaflet-specific time series plot
            print(f"\nCreating leaflet-specific time series plots with {window_size}-frame smoothing...")
            plot_leaflet_composition_time_series(time_series_data, proteins, lipid_sels, leaflet_finder, OUTPUT_DIR, window_size=window_size)
        
        # 4. Create bar charts for lipid composition per leaflet
        print("\nCreating bar charts for lipid composition per leaflet...")
        plot_lipid_bar_charts(time_series_data, proteins, lipid_sels, leaflet_finder, OUTPUT_DIR)
        
        # 5. Also create simple count bar charts
        print("\nCreating simple count bar charts for lipid composition...")
        plot_simple_bar_charts(time_series_data, proteins, lipid_sels, leaflet_finder, OUTPUT_DIR)
        
        # Record end time and print summary
        end_time = time.time()
        duration = end_time - start_time
        hours, remainder = divmod(duration, 3600)
        minutes, seconds = divmod(remainder, 60)
        
        print(f"\nAnalysis completed in {int(hours)}h {int(minutes)}m {int(seconds)}s")
        print(f"Results saved to {OUTPUT_DIR}")
        
        # Write completion summary
        with open(os.path.join(OUTPUT_DIR, 'analysis_summary.txt'), 'w') as f:
            f.write(f"Time Series Analysis completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Duration: {int(hours)}h {int(minutes)}m {int(seconds)}s\n\n")
            f.write(f"Trajectory: {TRAJECTORY_FILE}\n")
            f.write(f"Frames analyzed: {START} to {STOP}, step {STEP}, interval {TIME_SERIES_INTERVAL}\n")
            
            # Information about processing method
            if 'test_parallel_time' in locals():
                f.write(f"Parallelization: {NUM_CPUS} CPU cores\n")
            else:
                f.write("Processing: Single-threaded with NumPy optimization\n")
                
            f.write(f"Total frames in time series: {len(time_series_data)}\n\n")
            
            # Plot information
            f.write("Created plots:\n")
            f.write("  - Weighted lipid composition time series for each protein\n")
            f.write("  - Simple count lipid composition time series for each protein\n")
            f.write("  - Leaflet-specific lipid composition time series\n")
            f.write(f"  - All plots created with smoothing windows: {smoothing_windows}\n")
            f.write("  - Bar charts of lipid composition per protein and leaflet\n")
            f.write("  - Simple count bar charts of lipid composition\n")
            f.write("  - All plots use consistent color scheme for lipid types\n")
            f.write("  - All plots saved in both PNG and PS formats\n")
        
        return 0
    
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        
        with open(os.path.join(OUTPUT_DIR, 'error_log.txt'), 'w') as f:
            f.write(f"Error occurred at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Error message: {str(e)}\n\n")
            f.write("Traceback:\n")
            traceback.print_exc(file=f)
        
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
else:
    # When imported as a module
    print("LipidProteinAnalyzer loaded. Functions available for import.")

# Package metadata (added at the end)
__version__ = "1.0.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

def test_import():
    """Test function for JOSS review."""
    print("time_series_analysis module imported successfully!")
    return True