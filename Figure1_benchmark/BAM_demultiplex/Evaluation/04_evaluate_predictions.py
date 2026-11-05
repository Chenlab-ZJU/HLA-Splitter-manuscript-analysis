# Evaluate existing predictions only; no demultiplexing commands are executed.
# Legacy Accuracy/Accuracy_Singlets is correct / all input cells, not correct / assigned.
import os
import pandas as pd
import re

def normalise_barcodes(values):
    """Match 10x barcodes with or without a GEM suffix; reject collisions at joins."""
    return values.astype(str).str.split('-').str[0]


def locate_time_log(directory, *legacy_names):
    """Prefer the unified log while retaining support for existing run directories."""
    candidates = [os.path.join(directory, name) for name in ('time.log', *legacy_names)]
    return next((path for path in candidates if os.path.exists(path)), candidates[0])

def parse_time_log(log_file_path):
    """
    Parses the output of /usr/bin/time -v to extract relevant metrics.
    """
    time_metrics = {}
    if not os.path.exists(log_file_path):
        return time_metrics

    with open(log_file_path, 'r') as f:
        content = f.read()

    # Regex to extract key metrics
    # Wall clock time (elapsed time)
    elapsed = re.findall(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)', content)
    if elapsed:
        seconds = 0.0
        for value in elapsed:
            parts = [float(p) for p in value.split(':')]
            seconds += sum(part * 60 ** i for i, part in enumerate(reversed(parts)))
        time_metrics['Wall_Time_Sec'] = seconds
    
    # Maximum resident set size (peak memory usage)
    memory = re.findall(r'Maximum resident set size \(kbytes\): (\d+)', content)
    if memory:
        time_metrics['Peak_Memory_MB'] = max(map(int, memory)) / 1024

    # CPU usage percentage
    cpu = re.findall(r'Percent of CPU this job got: (\d+)%', content)
    if cpu:
        time_metrics['CPU_Percent'] = sum(map(int, cpu)) / len(cpu)

    return time_metrics

def map_clusters_to_ground_truth(tool_df, ground_truth_df, 
                                 pred_col='Predicted_Sample_ID', 
                                 true_col='True_Sample_ID', 
                                 barcode_col='Cell_Barcode'):
    """
    Based on the majority voting method, the random cluster names of the genotype-deconvolution tool are mapped to the actual sample IDs.
    """
    ignore_labels = ['doublet', 'unassigned']
    mapped_df = tool_df.copy()
    mapped_df[barcode_col] = normalise_barcodes(mapped_df[barcode_col])
    ground_truth_df = ground_truth_df.copy()
    ground_truth_df[barcode_col] = normalise_barcodes(ground_truth_df[barcode_col])
    
    mapped_df[pred_col] = mapped_df[pred_col].astype(str)
    singlets_mask = ~mapped_df[pred_col].isin(ignore_labels)
    singlets_only = mapped_df[singlets_mask].copy()
    
    temp_merge = pd.merge(singlets_only, ground_truth_df, on=barcode_col, how='inner', validate='one_to_one')
    
    # Group statistics and extract the majority to build a dictionary
    cluster_map = {}
    for cluster_id, group in temp_merge.groupby(pred_col):
        if not group.empty:
            majority_true_id = group[true_col].mode()[0]
            cluster_map[str(cluster_id)] = str(majority_true_id)
            
    original_predictions = mapped_df[pred_col].copy()
    mapped_df[pred_col] = mapped_df[pred_col].map(cluster_map)
    for label in ignore_labels:
        mapped_df.loc[original_predictions == label, pred_col] = label
        
    mapped_df[pred_col] = mapped_df[pred_col].fillna('unassigned')
    
    return mapped_df

def evaluate_demultiplexing(ground_truth_df, tool_output_df, tool_name):
    """
    Evaluates demultiplexing accuracy given ground truth and tool output.
    Assumes tool_output_df has 'Cell_Barcode' and 'Predicted_Sample_ID' (or similar).
    """
    ground_truth_df = ground_truth_df.copy()
    tool_output_df = tool_output_df.copy()
    for frame in (ground_truth_df, tool_output_df):
        frame['Cell_Barcode'] = normalise_barcodes(frame['Cell_Barcode'])
    merged_df = pd.merge(ground_truth_df, tool_output_df, on='Cell_Barcode', how='left', validate='one_to_one')
    total_ground_truth_cells = len(merged_df)
    
    # Standardize predicted sample IDs to string for comparison
    merged_df['Predicted_Sample_ID'] = merged_df['Predicted_Sample_ID'].astype(str)
    doublet_cells = merged_df[merged_df['Predicted_Sample_ID'].astype(str) == 'doublet']
    ambiguous_cells = merged_df[merged_df['Predicted_Sample_ID'].astype(str) == 'unassigned'] 
    assigned_singlets = merged_df[merged_df['Predicted_Sample_ID'].isin(ground_truth_df['True_Sample_ID'].unique())]
    
    correctly_assigned_singlets = assigned_singlets[
    assigned_singlets['True_Sample_ID'] == assigned_singlets['Predicted_Sample_ID']]

    num_total_cells = len(ground_truth_df)
    num_assigned_singlets = len(assigned_singlets)
    num_correctly_assigned_singlets = len(correctly_assigned_singlets)
    num_doublets = len(doublet_cells)
    num_ambiguous = len(ambiguous_cells)
    num_unassigned_by_tool = num_total_cells - num_assigned_singlets - num_doublets - num_ambiguous 

    # Preserve the original all-input denominator for the legacy metric.
    accuracy = num_correctly_assigned_singlets / num_total_cells if num_total_cells > 0 else 0

    # Per-sample metrics (Precision, Recall, F1-score)
    per_sample_metrics = {}
    for sample_id in ground_truth_df['True_Sample_ID'].unique():
        true_positives = correctly_assigned_singlets[correctly_assigned_singlets['True_Sample_ID'] == sample_id]
        true_sample_cells = ground_truth_df[ground_truth_df['True_Sample_ID'] == sample_id]
        predicted_sample_cells = assigned_singlets[assigned_singlets['Predicted_Sample_ID'] == sample_id]

        tp = len(true_positives)
        fp = len(predicted_sample_cells) - tp
        fn = len(true_sample_cells) - tp 

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

        per_sample_metrics[sample_id] = {
            'TP': tp, 'FP': fp, 'FN': fn,
            'Precision': precision, 'Recall': recall, 'F1_Score': f1_score
        }

    return {
        'Tool_Name': tool_name,
        'Total_Cells_Ground_Truth': num_total_cells,
        'Assigned_Singlets': num_assigned_singlets,
        'Correctly_Assigned_Singlets': num_correctly_assigned_singlets,
        'Accuracy_Singlets': accuracy,
        'Correct_Input_Fraction': accuracy,
        'Assigned_Singlet_Accuracy': num_correctly_assigned_singlets / num_assigned_singlets if num_assigned_singlets else None,
        'Cell_Recovery': num_assigned_singlets / num_total_cells if num_total_cells else 0,
        'Doublet_Calls': num_doublets,
        'Ambiguous_Calls': num_ambiguous,
        'Unassigned_By_Tool': num_unassigned_by_tool, 
        'Per_Sample_Metrics': per_sample_metrics
    }

# --- Main Evaluation Loop ---
if __name__ == "__main__":
    missing = [name for name in ('PROJECT_DIR', 'OUT_DIR', 'CELL_NUMS_STR') if not os.environ.get(name)]
    if missing:
        raise SystemExit('Source config.sh first; missing: ' + ', '.join(missing))
    PROJECT_DIR = os.environ.get("PROJECT_DIR")
    OUT_DIR = os.environ.get("OUT_DIR")
    # Split the space-separated pool names into a Python list.
    cell_nums = os.environ.get("CELL_NUMS_STR").split()
    
    benchmark_results = {}
    for cell_num in cell_nums:
        print(f"\n--- Evaluating results for simulated pool: {cell_num} ---")
        
        CURRENT_POOL_DIR = os.path.join(OUT_DIR, cell_num)
        GROUND_TRUTH_FILE = os.path.join(PROJECT_DIR, "True_sample_labels", f"{cell_num}.csv") 

        if not os.path.exists(GROUND_TRUTH_FILE):
            print(f"Error: Ground truth file {GROUND_TRUTH_FILE} not found. Skipping evaluation for {cell_num}.")
            continue
        
        # Standard .csv uses comma separator. If yours is tab-separated, change back to sep='\t'
        ground_truth_df = pd.read_csv(GROUND_TRUTH_FILE, sep=',')
        ground_truth_df = ground_truth_df[['barcode_select', 'sample']].rename(columns={'barcode_select': 'Cell_Barcode', 'sample': 'True_Sample_ID'})
        ground_truth_df['Cell_Barcode'] = normalise_barcodes(ground_truth_df['Cell_Barcode'])
        if ground_truth_df['Cell_Barcode'].duplicated().any():
            raise ValueError('Ground-truth barcodes collide after suffix normalization: ' + GROUND_TRUTH_FILE)
        print(f"Loaded {len(ground_truth_df)} cells from ground truth.")

        # --- Demuxlet Evaluation ---
        demuxlet_dir = os.path.join(CURRENT_POOL_DIR, "Demuxlet_res")
        demuxlet_output_file = os.path.join(demuxlet_dir, "demuxlet_results.best") 
        demuxlet_time_log = os.path.join(demuxlet_dir, "time.log")

        if os.path.exists(demuxlet_output_file):
            print("  Evaluating Demuxlet...")
            demuxlet_df = pd.read_csv(demuxlet_output_file, sep='\t')
            demuxlet_df = demuxlet_df[['BARCODE', 'DROPLET.TYPE', 'SNG.BEST.GUESS']].rename(columns={'BARCODE': 'Cell_Barcode', 'SNG.BEST.GUESS': 'Predicted_Sample_ID'})
            
            demuxlet_df['Predicted_Sample_ID'] = demuxlet_df['Predicted_Sample_ID'].astype(str)
            demuxlet_df.loc[demuxlet_df['DROPLET.TYPE'] == 'DBL', 'Predicted_Sample_ID'] = 'doublet'
            demuxlet_df.loc[demuxlet_df['DROPLET.TYPE'] == 'AMB', 'Predicted_Sample_ID'] = 'unassigned'

            demuxlet_metrics = evaluate_demultiplexing(ground_truth_df, demuxlet_df, "Demuxlet")
            demuxlet_metrics['Time_Metrics'] = parse_time_log(demuxlet_time_log)
            benchmark_results[f'{cell_num}_Demuxlet'] = demuxlet_metrics
            print(f"    Demuxlet Accuracy: {demuxlet_metrics['Accuracy_Singlets']:.2f}")
        else:
            print(f"  Demuxlet output file not found: {demuxlet_output_file}")


        # --- Vireo_res_with_genotypes Evaluation ---
        vireo_dir = os.path.join(CURRENT_POOL_DIR, "Vireo_res_with_genotypes")
        # Renamed variable to vireo_output_file to match checking logic
        vireo_output_file = os.path.join(vireo_dir, "vireo_results/donor_ids.tsv")
        if not os.path.exists(vireo_output_file):
            vireo_output_file = os.path.join(vireo_dir, "vireo_results_genotypes/donor_ids.tsv")
        vireo_time_log = locate_time_log(vireo_dir, 'vireo_run_time.log', 'vireo_time.log')

        if os.path.exists(vireo_output_file):
            print("  Evaluating Vireo (With Genotypes)...")
            vireo_df = pd.read_csv(vireo_output_file, sep='\t')
            label_col = os.environ.get('VIREO_GT_LABEL_COLUMN', 'best_singlet')
            if label_col not in ('best_singlet', 'donor_id'):
                raise ValueError('VIREO_GT_LABEL_COLUMN must be best_singlet or donor_id')
            vireo_df = vireo_df[['cell', label_col]].rename(columns={'cell': 'Cell_Barcode', label_col: 'Predicted_Sample_ID'})

            vireo_metrics = evaluate_demultiplexing(ground_truth_df, vireo_df, "Vireo_with_genotypes")
            vireo_metrics['Time_Metrics'] = parse_time_log(vireo_time_log)
            benchmark_results[f'{cell_num}_Vireo_with_genotypes'] = vireo_metrics
            print(f"    Vireo Accuracy: {vireo_metrics['Accuracy_Singlets']:.2f}")
        else:
            print(f"  Vireo output file not found: {vireo_output_file}")

        # --- Vireo_res_without_genotypes Evaluation ---
        vireo_dir = os.path.join(CURRENT_POOL_DIR, "Vireo_res_without_genotypes")
        # Renamed variable to vireo_output_file
        vireo_output_file = os.path.join(vireo_dir, "vireo_results/donor_ids.tsv") 
        vireo_time_log = locate_time_log(vireo_dir, 'vireo_run_time.log', 'vireo_time.log')

        if os.path.exists(vireo_output_file):
            print("  Evaluating Vireo (Without Genotypes)...")
            vireo_df = pd.read_csv(vireo_output_file, sep='\t')
            vireo_df = vireo_df[['cell', 'donor_id']].rename(columns={'cell': 'Cell_Barcode', 'donor_id': 'Predicted_Sample_ID'})

            # Based on the majority voting method, the random cluster names of the genotype-deconvolution tool are mapped to the actual sample IDs.
            vireo_df = map_clusters_to_ground_truth(vireo_df, ground_truth_df, 
                pred_col='Predicted_Sample_ID', true_col='True_Sample_ID', barcode_col='Cell_Barcode')

            vireo_metrics = evaluate_demultiplexing(ground_truth_df, vireo_df, "Vireo_without_genotypes")
            vireo_metrics['Time_Metrics'] = parse_time_log(vireo_time_log)
            benchmark_results[f'{cell_num}_Vireo_without_genotypes'] = vireo_metrics
            print(f"  Vireo Accuracy: {vireo_metrics['Accuracy_Singlets']:.2f}")
        else:
            print(f"  Vireo output file not found: {vireo_output_file}")

        # --- Souporcell Evaluation ---
        souporcell_dir = os.path.join(CURRENT_POOL_DIR, "Souporcell_res")
        souporcell_output_file = os.path.join(souporcell_dir, "clusters.tsv") 
        souporcell_time_log = locate_time_log(souporcell_dir, 'souporcell_time.log')

        if os.path.exists(souporcell_output_file):
            print("  Evaluating Souporcell...")
            souporcell_df = pd.read_csv(souporcell_output_file, sep='\t')
            souporcell_df = souporcell_df[['barcode', 'status', 'assignment']].rename(columns={'barcode': 'Cell_Barcode', 'assignment': 'Predicted_Sample_ID'})
            
            souporcell_df['Predicted_Sample_ID'] = souporcell_df['Predicted_Sample_ID'].astype(str)
            souporcell_df.loc[souporcell_df['status'] == 'doublet', 'Predicted_Sample_ID'] = 'doublet'
            souporcell_df.loc[souporcell_df['status'] == 'unassigned', 'Predicted_Sample_ID'] = 'unassigned'

            # Based on the majority voting method, the random cluster names of the genotype-deconvolution tool are mapped to the actual sample IDs.
            souporcell_df = map_clusters_to_ground_truth(souporcell_df, ground_truth_df, 
                pred_col='Predicted_Sample_ID', true_col='True_Sample_ID', barcode_col='Cell_Barcode')

            souporcell_metrics = evaluate_demultiplexing(ground_truth_df, souporcell_df, "Souporcell")
            souporcell_metrics['Time_Metrics'] = parse_time_log(souporcell_time_log)
            benchmark_results[f'{cell_num}_Souporcell'] = souporcell_metrics
            print(f"    Souporcell Accuracy: {souporcell_metrics['Accuracy_Singlets']:.2f}")
        else:
            print(f"  Souporcell output file not found: {souporcell_output_file}")


        # --- HLA-Splitter Evaluation ---
        hla_splitter_dir = os.path.join(CURRENT_POOL_DIR, "HLA_Splitter_res")
        hla_splitter_output_file = os.path.join(hla_splitter_dir, "HLA_demultiplex.csv") 
        hla_splitter_time_log = locate_time_log(hla_splitter_dir, 'HLA_Splitter_time.log')

        if os.path.exists(hla_splitter_output_file):
            print("  Evaluating HLA-Splitter...")
            # Standard .csv uses comma separator
            hla_splitter_df = pd.read_csv(hla_splitter_output_file, sep=',')
            label_col = 'Sample_prediction' if 'Sample_prediction' in hla_splitter_df else 'Prediction'
            hla_splitter_df = hla_splitter_df[['Cell barcode', label_col]].rename(columns={'Cell barcode': 'Cell_Barcode', label_col: 'Predicted_Sample_ID'})

            hla_splitter_df['Predicted_Sample_ID'] = hla_splitter_df['Predicted_Sample_ID'].astype(str)
            hla_splitter_df['Cell_Barcode'] = normalise_barcodes(hla_splitter_df['Cell_Barcode'])
            hla_splitter_df.loc[hla_splitter_df['Predicted_Sample_ID'] == 'Doublet', 'Predicted_Sample_ID'] = 'doublet'
            hla_splitter_df.loc[hla_splitter_df['Predicted_Sample_ID'] == 'Undefined', 'Predicted_Sample_ID'] = 'unassigned'

            hla_splitter_metrics = evaluate_demultiplexing(ground_truth_df, hla_splitter_df, "HLA-Splitter")
            hla_splitter_metrics['Time_Metrics'] = parse_time_log(hla_splitter_time_log)
            benchmark_results[f'{cell_num}_HLA-Splitter'] = hla_splitter_metrics
            print(f"    HLA-Splitter Accuracy: {hla_splitter_metrics['Accuracy_Singlets']:.2f}")
        else:
            print(f"  HLA-Splitter output file not found: {hla_splitter_output_file}")

    if not benchmark_results:
        raise SystemExit('No matching predictions found; existing summaries were not modified.')

    # --- Print summary results ---
    print("\n--- Benchmark Summary ---")
    csv_data = []

    for key, metrics in benchmark_results.items():
        print(f"\nResults for {key}:")
        print(f"  Total Cells: {metrics['Total_Cells_Ground_Truth']}")
        print(f"  Assigned Singlets: {metrics['Assigned_Singlets']}")
        print(f"  Correctly Assigned Singlets: {metrics['Correctly_Assigned_Singlets']}")
        print(f"  Accuracy (overall): {metrics['Accuracy_Singlets']:.4f}")
        print(f"  Doublet Calls: {metrics['Doublet_Calls']}")
        print(f"  Ambiguous Calls: {metrics['Ambiguous_Calls']}")
        print(f"  Unassigned by Tool: {metrics['Unassigned_By_Tool']}")
        
        wall_time = 'N/A'
        peak_mem = 'N/A'
        cpu_pct = 'N/A'

        if 'Time_Metrics' in metrics and metrics['Time_Metrics']:
            wall_time = metrics['Time_Metrics'].get('Wall_Time_Sec', 'N/A')
            peak_mem = metrics['Time_Metrics'].get('Peak_Memory_MB', 'N/A')
            cpu_pct = metrics['Time_Metrics'].get('CPU_Percent', 'N/A')
            
            print(f"  Wall Time: {wall_time} seconds")
            if peak_mem != 'N/A':
                print(f"  Peak Memory: {peak_mem:.2f} MB")
            else:
                print(f"  Peak Memory: {peak_mem} MB")
            print(f"  CPU Percent: {cpu_pct}%")
        else:
            print("  Time metrics not available.")

        row_data = {
            'Run_ID': key,
            'Tool_Name': metrics.get('Tool_Name', 'N/A'),
            'Total_Cells': metrics.get('Total_Cells_Ground_Truth', 'N/A'),
            'Accuracy': metrics.get('Accuracy_Singlets', 'N/A'),
            'Correct_Input_Fraction': metrics['Correct_Input_Fraction'],
            'Assigned_Singlet_Accuracy': metrics['Assigned_Singlet_Accuracy'],
            'Cell_Recovery': metrics['Cell_Recovery'],
            'Vireo_Genotype_Label_Policy': os.environ.get('VIREO_GT_LABEL_COLUMN', 'best_singlet'),
            'Assigned_Singlets': metrics.get('Assigned_Singlets', 'N/A'),
            'Correctly_Assigned': metrics.get('Correctly_Assigned_Singlets', 'N/A'),
            'Doublets': metrics.get('Doublet_Calls', 'N/A'),
            'Ambiguous': metrics.get('Ambiguous_Calls', 'N/A'),
            'Unassigned': metrics.get('Unassigned_By_Tool', 'N/A'),
            'Wall_Time_Sec': wall_time,
            'Peak_Memory_MB': peak_mem,
            'CPU_Percent': cpu_pct
        }
        csv_data.append(row_data)

    # 1. save as CSV file (flattened structure for easy analysis)
    summary_df = pd.DataFrame(csv_data)
    output_csv_path = os.path.join(PROJECT_DIR, "benchmark_summary.csv")
    summary_df.to_csv(output_csv_path, index=False)
    print(f"\nDetailed benchmark results saved to CSV: {output_csv_path}")

    # 2. save as JSON file (preserve deep nested data structure)
    output_summary_path = os.path.join(PROJECT_DIR, "benchmark_summary.json")
    import json
    with open(output_summary_path, 'w') as f:
        json.dump(benchmark_results, f, indent=4)
    print(f"Detailed benchmark results saved to JSON: {output_summary_path}")
