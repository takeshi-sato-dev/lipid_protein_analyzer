#!/usr/bin/env python
"""
Generate test data for LipidProteinAnalyzer tests.
Creates minimal but realistic topology and trajectory files.
"""

import os
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def create_test_universe(n_proteins=2, n_lipids_per_type=10, n_frames=10):
    """
    Create a test universe with proteins and lipids.
    
    Parameters
    ----------
    n_proteins : int
        Number of proteins to create
    n_lipids_per_type : int
        Number of lipids per lipid type
    n_frames : int
        Number of trajectory frames
    
    Returns
    -------
    universe : MDAnalysis.Universe
        Test universe with topology
    trajectory_data : list
        List of position arrays for each frame
    """
    
    # Lipid types to include
    lipid_types = ['POPC', 'POPE', 'POPS', 'CHOL']
    n_lipid_types = len(lipid_types)
    n_total_lipids = n_lipids_per_type * n_lipid_types
    
    # Atoms per molecule (Martini coarse-grained)
    atoms_per_protein = 10  # Simplified protein
    atoms_per_lipid = 5     # Martini lipid (NC3, PO4, GL1, GL2, C1A)
    
    # Calculate total atoms
    n_protein_atoms = n_proteins * atoms_per_protein
    n_lipid_atoms = n_total_lipids * atoms_per_lipid
    n_atoms = n_protein_atoms + n_lipid_atoms
    
    print(f"Creating test system:")
    print(f"  - {n_proteins} proteins ({n_protein_atoms} atoms)")
    print(f"  - {n_total_lipids} lipids ({n_lipid_atoms} atoms)")
    print(f"  - Total: {n_atoms} atoms")
    
    # Build topology arrays
    atom_names = []
    atom_types = []
    resnames = []
    resids = []
    segids = []
    
    current_resid = 0
    
    # Add proteins
    for p in range(n_proteins):
        current_resid += 1
        for a in range(atoms_per_protein):
            if a == 0:
                atom_names.append('BB')  # Backbone bead
                atom_types.append('P1')
            else:
                atom_names.append(f'SC{a}')  # Sidechain
                atom_types.append('P2')
            
            resnames.append('PRO')
            resids.append(current_resid)
            segids.append(f'PRO{p}')
    
    # Add lipids
    for lipid_type in lipid_types:
        for l in range(n_lipids_per_type):
            current_resid += 1
            
            # Martini lipid beads
            if lipid_type == 'CHOL':
                beads = [('ROH', 'SP1'), ('R1', 'SC1'), ('R2', 'SC1'), 
                        ('R3', 'SC1'), ('R4', 'SC1')]
            else:
                # Phospholipid
                if lipid_type == 'POPC':
                    head = ('NC3', 'Q0')  # Choline
                elif lipid_type == 'POPE':
                    head = ('NH3', 'Qd')  # Ethanolamine
                elif lipid_type == 'POPS':
                    head = ('CNO', 'Qa')  # Serine
                else:
                    head = ('NC3', 'Q0')  # Default
                
                beads = [
                    head,
                    ('PO4', 'Qa'),  # Phosphate
                    ('GL1', 'Na'),  # Glycerol
                    ('GL2', 'Na'),  # Glycerol
                    ('C1A', 'C1'),  # Tail
                ]
            
            for bead_name, bead_type in beads:
                atom_names.append(bead_name)
                atom_types.append(bead_type)
                resnames.append(lipid_type)
                resids.append(current_resid)
                segids.append('MEMB')
    
    # Create universe
    u = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=current_resid,
        atom_resindex=np.array([r-1 for r in resids]),
        trajectory=True
    )
    
    # Add topology attributes
    u.add_TopologyAttr('name', atom_names)
    u.add_TopologyAttr('type', atom_types)
    u.add_TopologyAttr('resname', resnames)
    u.add_TopologyAttr('resid', resids)
    u.add_TopologyAttr('segid', segids)
    
    # Generate trajectory frames
    box_size = 100.0  # Angstroms
    trajectory_data = []
    
    for frame_idx in range(n_frames):
        positions = np.zeros((n_atoms, 3))
        atom_idx = 0
        
        # Position proteins (membrane-embedded)
        for p in range(n_proteins):
            # Distribute proteins in XY plane
            protein_x = 30 + p * 40
            protein_y = 50
            protein_z = 50  # Middle of membrane
            
            for a in range(atoms_per_protein):
                # Arrange in a helix-like structure
                angle = a * np.pi / 5
                x = protein_x + 5 * np.cos(angle)
                y = protein_y + 5 * np.sin(angle)
                z = protein_z + a * 2
                
                positions[atom_idx] = [x, y, z]
                atom_idx += 1
        
        # Position lipids in two leaflets
        lipids_per_leaflet = n_total_lipids // 2
        
        for leaflet in range(2):
            z_base = 70 if leaflet == 0 else 30  # Upper and lower leaflets
            
            for lipid_idx in range(lipids_per_leaflet):
                # Grid arrangement in XY plane
                grid_size = int(np.sqrt(lipids_per_leaflet))
                x_idx = lipid_idx % grid_size
                y_idx = lipid_idx // grid_size
                
                x_base = 10 + x_idx * (box_size - 20) / grid_size
                y_base = 10 + y_idx * (box_size - 20) / grid_size
                
                # Position lipid beads
                for b in range(atoms_per_lipid):
                    x = x_base + np.random.randn() * 1
                    y = y_base + np.random.randn() * 1
                    z = z_base + b * 2  # Stack beads vertically
                    
                    positions[atom_idx] = [x, y, z]
                    atom_idx += 1
        
        # Add random motion for dynamics
        if frame_idx > 0:
            positions += np.random.randn(n_atoms, 3) * 0.5
        
        trajectory_data.append(positions)
    
    # Set initial positions
    u.atoms.positions = trajectory_data[0]
    u.dimensions = [box_size, box_size, box_size, 90, 90, 90]
    
    return u, trajectory_data


def write_test_files(universe, trajectory_data, output_dir='test_data'):
    """
    Write test topology and trajectory files.
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe with topology
    trajectory_data : list
        List of position arrays for each frame
    output_dir : str
        Output directory for files
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Write PSF file
    psf_file = os.path.join(output_dir, 'test_system.psf')
    universe.atoms.write(psf_file)
    print(f"✓ Created topology: {psf_file}")
    
    # Write XTC trajectory
    xtc_file = os.path.join(output_dir, 'test_trajectory.xtc')
    n_atoms = len(universe.atoms)
    
    with XTCWriter(xtc_file, n_atoms=n_atoms) as writer:
        for frame_idx, positions in enumerate(trajectory_data):
            universe.atoms.positions = positions
            universe.trajectory.ts.frame = frame_idx
            universe.trajectory.ts.time = frame_idx * 1000  # ps
            writer.write(universe.atoms)
    
    print(f"✓ Created trajectory: {xtc_file}")
    
    # Write info file
    info_file = os.path.join(output_dir, 'test_info.txt')
    with open(info_file, 'w') as f:
        f.write("Test Data Information\n")
        f.write("====================\n\n")
        f.write(f"Topology: {psf_file}\n")
        f.write(f"Trajectory: {xtc_file}\n")
        f.write(f"Frames: {len(trajectory_data)}\n")
        f.write(f"Atoms: {n_atoms}\n")
        f.write(f"Box: {universe.dimensions[:3]} Å\n")
        f.write("\nResidues:\n")
        
        # Count residues by type
        resname_counts = {}
        for resname in universe.atoms.resnames:
            resname_counts[resname] = resname_counts.get(resname, 0) + 1
        
        for resname, count in resname_counts.items():
            f.write(f"  {resname}: {count} atoms\n")
    
    print(f"✓ Created info file: {info_file}")
    
    return psf_file, xtc_file


def verify_test_data(psf_file, xtc_file):
    """
    Verify that test data is valid and usable.
    
    Parameters
    ----------
    psf_file : str
        Path to PSF file
    xtc_file : str
        Path to XTC file
    
    Returns
    -------
    bool
        True if data is valid
    """
    
    print("\nVerifying test data...")
    
    try:
        # Load universe
        u = mda.Universe(psf_file, xtc_file)
        print(f"✓ Loaded universe: {len(u.atoms)} atoms, {len(u.trajectory)} frames")
        
        # Check proteins
        proteins = u.select_atoms("resname PRO")
        print(f"✓ Found proteins: {len(proteins)} atoms")
        
        # Check lipids
        lipids = u.select_atoms("resname POPC POPE POPS CHOL")
        print(f"✓ Found lipids: {len(lipids)} atoms")
        
        # Check lipid types
        for lipid_type in ['POPC', 'POPE', 'POPS', 'CHOL']:
            selection = u.select_atoms(f"resname {lipid_type}")
            if len(selection) > 0:
                print(f"  - {lipid_type}: {len(selection)} atoms")
        
        # Check trajectory
        print(f"✓ Trajectory: {len(u.trajectory)} frames")
        
        # Check leaflet separation
        po4_atoms = u.select_atoms("name PO4")
        if len(po4_atoms) > 0:
            z_coords = po4_atoms.positions[:, 2]
            if len(np.unique(z_coords > 50)) == 2:
                print("✓ Leaflets separable (lipids in upper/lower regions)")
        
        # Import and test analysis functions
        try:
            import time_series_analysis as tsa
            
            # Override file paths
            original_psf = tsa.TOPOLOGY_FILE
            original_xtc = tsa.TRAJECTORY_FILE
            
            tsa.TOPOLOGY_FILE = psf_file
            tsa.TRAJECTORY_FILE = xtc_file
            
            # Try loading
            test_u = tsa.load_universe()
            print("✓ Compatible with time_series_analysis.load_universe()")
            
            # Try leaflet identification
            leaflet0, leaflet1, L = tsa.identify_lipid_leaflets(test_u)
            print("✓ Leaflet identification works")
            
            # Restore original paths
            tsa.TOPOLOGY_FILE = original_psf
            tsa.TRAJECTORY_FILE = original_xtc
            
        except Exception as e:
            print(f"⚠ Warning: Could not test with analysis module: {e}")
        
        return True
        
    except Exception as e:
        print(f"✗ Verification failed: {e}")
        return False


def main():
    """Main function to generate test data."""
    
    print("="*60)
    print("LipidProteinAnalyzer - Test Data Generator")
    print("="*60)
    
    # Parameters for test data
    n_proteins = 2
    n_lipids_per_type = 10  # 10 of each: POPC, POPE, POPS, CHOL
    n_frames = 10
    
    # Create test universe
    universe, trajectory_data = create_test_universe(
        n_proteins=n_proteins,
        n_lipids_per_type=n_lipids_per_type,
        n_frames=n_frames
    )
    
    # Write files
    psf_file, xtc_file = write_test_files(universe, trajectory_data)
    
    # Verify
    if verify_test_data(psf_file, xtc_file):
        print("\n✅ Test data generated successfully!")
        print("\nYou can now run:")
        print("  pytest tests/")
        print("  python test_quick.py")
    else:
        print("\n❌ Test data generation failed!")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())