import uproot
import pandas as pd
import numpy as np
import os
import argparse
import sys
from tqdm import tqdm

def load_calibration_db(calib_root_path):
    """
    Loads the calibration database into memory.
    Returns:
        - id_df: DataFrame containing the ID columns and an 'index' column.
        - plane_array: Numpy array of the 'plane' branch (N_events, 5).
        - prob_array: Numpy array of the 'prob' branch (N_events, 5).
    """
    print(f"--- Loading calibration database from: {calib_root_path} ---")
    print("This may take a moment depending on file size...")

    try:
        with uproot.open(calib_root_path) as f:
            tree = f["calibrated"]
            
            id_branches = ["run", "subrun", "gate", "physEvt"]
            id_df = tree.arrays(id_branches, library="pd")
            
            id_df = id_df.rename(columns={
                "run": "ev_run",
                "subrun": "ev_subrun",
                "gate": "ev_gate",
                "physEvt": "physEvtNum"
            })
            
            
            for col in id_df.columns:
                id_df[col] = id_df[col].astype(np.int64) # Ensure types are consistent (int64) to avoid merge mismatches

            id_df["db_idx"] = np.arange(len(id_df))
            
            payload = tree.arrays(["planes", "probs"], library="np")
            plane_array = payload["planes"]
            prob_array = payload["probs"]
            
            print(f"Database loaded. Total calibrated events: {len(id_df)}")
            return id_df, plane_array, prob_array

    except Exception as e:
        print(f"Error loading calibration file: {e}")
        sys.exit(1)

def process_file(original_path, output_dir, calib_id_df, calib_planes, calib_probs):
    """
    Process a single ROOT file: matches events and writes the friend file.
    """
    filename = os.path.basename(original_path)
    
    base_name = os.path.splitext(filename)[0]
    output_filename = f"{base_name}_calibrated_probs.root"
    output_path = os.path.join(output_dir, output_filename)

    try:
        with uproot.open(original_path) as f:
            if "MasterAnaDev" not in f:
                print(f"Skipping {filename}: 'MasterAnaDev' tree not found.")
                return

            tree = f["MasterAnaDev"]
            n_entries = tree.num_entries
            
            if n_entries == 0:
                print(f"Skipping {filename}: Empty tree.")
                return

            target_ids = tree.arrays(["ev_run", "ev_subrun", "ev_gate", "physEvtNum"], library="pd")
            
            for col in target_ids.columns:
                target_ids[col] = target_ids[col].astype(np.int64)

            merged = pd.merge(target_ids, calib_id_df, on=["ev_run", "ev_subrun", "ev_gate", "physEvtNum"], how="left") #merging the Original IDs (Left) with the Calib DB IDs (Right) preserving the order and length of the Original IDs.
            
            out_planes = np.full((n_entries, 5), -1, dtype=np.int32)
            out_probs = np.full((n_entries, 5), -1.0, dtype=np.float32)
            
            # Find which rows had a match
            # 'db_idx' will be NaN for non-matched rows
            match_mask = merged["db_idx"].notna()
            
            if match_mask.any():
                # Get the destination indices (in the friend file)
                dest_indices = np.where(match_mask)[0]
                
                # Get the source indices (from the calibration DB)
                src_indices = merged.loc[match_mask, "db_idx"].astype(int).values
                
                # Copy data from DB to Output
                out_planes[dest_indices] = calib_planes[src_indices]
                out_probs[dest_indices] = calib_probs[src_indices]
            
            output_data = {
                "plane": out_planes,
                "prob": out_probs
            }
            
            with uproot.recreate(output_path) as out_f:
                out_f.mktree('calibrated', output_data)

    except Exception as e:
        print(f"Error processing {filename}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Generate Calibrated Friend ROOT Files")
    parser.add_argument("-i", "--input-list", required=True, help="Path to .txt file containing list of original ROOT files")
    parser.add_argument("-c", "--calib-file", required=True, help="Path to large ROOT file containing calibrated outputs")
    parser.add_argument("-o", "--output-dir", required=True, help="Directory to save the new friend files")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        print(f"Creating output directory: {args.output_dir}")
        os.makedirs(args.output_dir)
        
    id_df, plane_arr, prob_arr = load_calibration_db(args.calib_file)
    
    with open(args.input_list, 'r') as f:
        file_list = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(file_list)} files to process.")
    
    for file_path in tqdm(file_list, desc="Processing Files", unit="file"):
        process_file(file_path, args.output_dir, id_df, plane_arr, prob_arr)
        
    print("Done!")

if __name__ == "__main__":
    main()