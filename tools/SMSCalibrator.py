#This script trains a structures matrix scaling (or other specified type) based calibrator for MINERvA DANN predictions
#Based on https://arxiv.org/pdf/2511.03685

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import log_loss
from tqdm import tqdm
import os, joblib, argparse

import sys
from probmetrics.calibrators import get_calibrator, SMSCalibrator, Calibrator
import dirichletcal
from typing import Optional, List


OUTPUT_DIR = "calibration_results"

# Col 5 (index 4) = True Label
# Col 6-179 (indices 5-179) = Probabilities
LABEL_COL_IDX = 4 #Col 5 of the DANN raw prediction csv-like text file (index 4) corresponds to the true plane
PROB_COL_START = 6 #Col 7 (index 6) is the start of the actual plane probabilites, with this being plane 0
PROB_COL_END = 220 #Col 221 (index 220) is the end of the actual plane probabilities corresponding to plane 214

#Note, in some earlier iterations of the DANN there were only 174 planes (No DSECAL), if using these for whatever reason make sure to change PROB_COL_END appropriately

#Creating output directory
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def load_data_buffered(filepath, subsample_rate=1.0):
    """
    Because the raw DANN prediction files is HUGE, I need to use a buffered and chunked approach to prevent this from crashing due to RAM overuse on file loading
    Additionally I apply a cut here, I remove events that don't have an uncalibrated ML prediction since there's no point incuding them and contaminating the data we use to train our calibrator
    Also performs random subsampling
    """
    print(f"Reading large file: {filepath}")
    print(f"Subsampling rate: {subsample_rate * 100}%")
    
    labels_list = []
    probs_list = []
    
    file_size = os.path.getsize(filepath)
    with open(filepath, 'r') as f:
        with tqdm(total=file_size, unit='B', unit_scale=True, desc="Loading Data") as pbar: #Using tqdm to have some idea of progress as this can take a while
            
            #Using chunksize to break up input file reading
            reader = pd.read_csv(f, header=None, chunksize=100000)  # Process 100k rows at a time
            
            for chunk in reader:
                raw_labels = chunk.iloc[:, LABEL_COL_IDX].values
                cnn_labels = chunk.iloc[:, LABEL_COL_IDX+1].values
                valid_mask = (cnn_labels > 0) # A mask to keep only the entries that have a CNN prediction - the CNN prediction is stored the in row after the truth plane

                
                # Using another mask to subsample to save RAM
                if subsample_rate < 1.0:
                    # Generate a random mask of the same size
                    # Works by creating an array of numbers randomly distrubuted between 0-1 equal in length to the number of text lines/entries in this chunk of data and then evalues each of these array values with < subsample_rate 
                    random_mask = np.random.rand(len(chunk)) < subsample_rate

                    final_mask = valid_mask & random_mask # Combined mask
                else: #If not subsampling can skip the np.random stuff
                    final_mask = valid_mask

                labels = raw_labels[final_mask].astype(np.int64) # Apply filter masks
                probs = chunk.iloc[final_mask, PROB_COL_START : PROB_COL_END + 1].values.astype(np.float32)
                
                if len(labels) > 0:
                    labels_list.append(labels)
                    probs_list.append(probs)
                
                pbar.update(f.tell() - pbar.n)
                #break

    print("Concatenating chunks...")
    if not labels_list:
        raise ValueError("No data found! Check your filter, subsample rate or input file.")
        
    y = np.concatenate(labels_list)
    X_probs = np.concatenate(probs_list)
    
    # Force garbage collection
    del labels_list
    del probs_list
    import gc
    gc.collect()
    
    print(f"Data Loaded. Shape: {X_probs.shape} (Subsampled)")
    return X_probs, y

# ==========================================
#This ended up not being needed, the probmetrics package actually does this implicitly
# ==========================================
def probs_to_logits(probs, epsilon=1e-15):
    """
    Inverts Softmax to retrieve 'Pseudo-Logits'.
    x = ln(p)
    """
    print("Converting probabilities to pseudo-logits...")
    probs_clipped = np.clip(probs, epsilon, 1.0) # Clipping probabilities to avoid ln(0)
    logits = np.log(probs_clipped)
    return logits

def plot_multirank_reliability(probs, y_true, title, filename, n_ranks=3, n_bins=10):
    """
    Plots N reliability diagrams for the top N predicted ranks.
    """
    print(f"Generating reliability diagram: {title}")
    
    plt.figure(figsize=(6 * n_ranks, 5))
    
    # Get the indices of the top N predictions
    # argsort sorts ascending, so we take the last n_ranks columns and reverse
    top_preds_indices = np.argsort(probs, axis=1)[:, -n_ranks:][:, ::-1]
    for rank in range(n_ranks):
        plt.subplot(1, n_ranks, rank + 1)
        
        # Get the class index predicted at this rank
        rank_class_indices = top_preds_indices[:, rank]
        
        # Get the probability assigned to this rank
        # Use numpy advanced indexing to grab the prob for the specific class at this rank
        confidences = probs[np.arange(len(probs)), rank_class_indices]

        # Determine if this rank was the correct class
        accuracies = (rank_class_indices == y_true).astype(float)
        
        # --- Binning ---
        bin_boundaries = np.linspace(0, 1, n_bins + 1)
        bin_lowers = bin_boundaries[:-1]
        bin_uppers = bin_boundaries[1:]
        
        ece = 0.0
        bin_accs = []
        bin_confs = []
        
        for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
            in_bin = (confidences > bin_lower) & (confidences <= bin_upper)
            prop_in_bin = in_bin.mean()
            
            if prop_in_bin > 0:
                accuracy = accuracies[in_bin].mean()
                avg_conf = confidences[in_bin].mean()
                
                ece += np.abs(avg_conf - accuracy) * prop_in_bin
                
                bin_accs.append(accuracy)
                bin_confs.append(avg_conf)
            else:
                bin_accs.append(np.nan)
                bin_confs.append(np.nan)

        # Plot - reproducing reliability plot style often seen in ML papers
        plt.plot([0, 1], [0, 1], "k--", label="Perfectly Calibrated")
        plt.scatter(bin_confs, bin_accs, color='blue')
        plt.plot(bin_confs, bin_accs, color='blue', alpha=0.5)
        plt.bar(bin_lowers, bin_accs, width=1.0/n_bins, align='edge', color='blue', alpha=0.2, edgecolor='black')
        
        plt.xlabel("Confidence")
        plt.ylabel("Accuracy")
        nth = rank+1
        suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(nth if nth < 20 else nth % 10, 'th')
        plt.title(f"{nth}{suffix} Most Proable Plane (ECE: {ece:.4f})")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.ylim(0, 1/(rank+1))
        plt.xlim(0, 1/(rank+1))

    plt.suptitle(title, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(OUTPUT_DIR, filename))
    plt.close()

    print("Finished Plot Generation")

def plot_nth_plane_distribution(probs, n=1, num_classes=None, ax=None, label=None, color=None, save_filename=None, title=None):
    """
    Plots a histogram of the plane numbers that correspond to the n-th highest 
    probability for each event and saves it to a file.

    Parameters:
    - probs: (numpy.ndarray) Matrix of shape [N_events, N_classes] containing probabilities.
    - n: (int) The rank to plot. 1 = Most Probable, 2 = 2nd Most Probable, etc.
    - num_classes: (int) Total number of planes. If None, inferred from probs.
    - ax: (matplotlib.axes.Axes) Optional axis to plot on. If None, creates a new figure.
    - label: (str) Label for the legend.
    - color: (str) Color of the histogram line.
    - save_filename: (str) Full path/filename to save the plot (e.g., 'rank1_calibrated.pdf').
                     If None, the plot is not saved automatically.
    """
    
    if num_classes is None:
        num_classes = probs.shape[1]

    # Sort and pick the n-th highest index (-1 is highest, -2 is 2nd highest, etc.)
    nth_indices = np.argsort(probs, axis=1)[:, -n]

    # Create figure if axis not provided otherwise use existing figure, this allows drawing multiple plots overlaid, if ever desired
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.arange(num_classes + 1) - 0.5 # Center bins on integers
    
    ax.hist(nth_indices, bins=bins, histtype='step', linewidth=2, label=label, color=color)
    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n if n < 20 else n % 10, 'th')
    ax.set_xlabel('Plane Number')
    ax.set_ylabel(f'Number of Events ({n}{suffix} Choice)')
    if title is None:
        ax.set_title(f'Distribution of the {n}{suffix} Most Probable Plane')
    else:
        ax.set_title(title)
    ax.set_xlim(0, num_classes)
    
    if label:
        ax.legend()

    # Save file
    if save_filename:
        print(f"Saving plot to {save_filename}...")
        ax.figure.savefig(save_filename, bbox_inches='tight', dpi=300)
        
    return ax

def plot_prob_distributions(probs, n=1, num_bins=10, ax=None, label=None, color=None, save_filename=None, title=None):
    """
    Plots a histogram of the plane numbers that correspond to the n-th highest 
    probability for each event and optionally saves it to a file.

    Parameters:
    - probs: (numpy.ndarray) Matrix of shape [N_events, N_classes] containing probabilities.
    - n: (int) The probability to plot. 1 = Most Probable, 2 = 2nd Most Probable, etc.
    - num_classes: (int) Total number of planes. If None, inferred from probs.
    - ax: (matplotlib.axes.Axes) Optional axis to plot on. If None, creates a new figure.
    - label: (str) Label for the legend.
    - color: (str) Color of the histogram line.
    - save_filename: (str) Full path/filename to save the plot (e.g., 'rank1_calibrated.pdf').
                     If None, the plot is not saved automatically.
    """
    
    # Sort and pick the n-th highest probability (-1 is highest, -2 is 2nd highest, etc.)
    nth_prob = np.sort(probs, axis=1)[:, -n]

    # Create figure if axis not provided otherwise use existing figure, this allows drawing multiple plots overlaid, if ever desired
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.linspace(0, 1, num_bins + 1)
    num_evts = len(nth_prob)
    weights = np.full(num_evts, 1/num_evts)
    ax.hist(nth_prob, bins=bins, histtype='step', linewidth=2, label=label, color=color, weights=weights)
    ax.set_xlabel('CNN Probability Score')
    ax.set_ylabel(f'Number of Events')
    if title is None:
        ax.set_title('Event Probability Distribtuion')
    else:
        ax.set_title(title)
    
    if label:
        ax.legend()

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        ax.figure.savefig(save_filename, bbox_inches='tight', dpi=300)
        
    return ax

def plot_prob_distributions_data_vs_mc(mc_uncalib, mc_calib, data_uncalib, data_calib, n=1, num_bins=10, ax=None, label=None, save_filename=None, title=None):
    
    fig, (ax_main, ax_ratio, ax_ratio2) = plt.subplots(
        nrows=3, 
        sharex=True, 
        gridspec_kw={'height_ratios': [3, 1, 1]}, # Top plot is 3x taller
        figsize=(8, 8)
    )
    plt.xlim(0, 1/n)
    plt.subplots_adjust(hspace=0)

    plot_prob_distributions(mc_uncalib, n=n, num_bins=num_bins, ax=ax_main, color='#d2c271', label="MC Uncalibrated")
    plot_prob_distributions(mc_calib, n=n, num_bins=num_bins, ax=ax_main, color='#CA5454', label="MC Calibrated")
    plot_prob_distributions(data_uncalib, n=n, num_bins=num_bins, ax=ax_main, color='#40B0A6', label="Data Uncalibrated")
    plot_prob_distributions(data_calib, n=n, num_bins=num_bins, ax=ax_main, color='#9285d5', label="Data Calibrated")

    if n==1:
        ax_main.axvline(0.4, color='gray', linestyle='--', linewidth=1)
    
    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n if n < 20 else n % 10, 'th')
    if title is None:
        ax_main.set_title(f'Probability Distribution of the {n}{suffix} Most Probable Plane')

    ax_ratio.axhline(1, color='gray', linestyle='--')
    ax_ratio.set_ylim(0, 2)
    bins = np.linspace(0, 1, num_bins + 1)

    nth_prob_data_uncalib = np.sort(data_uncalib, axis=1)[:, -n]
    nth_prob_mc_uncalib = np.sort(mc_uncalib, axis=1)[:, -n]
    nth_prob_data_calib = np.sort(data_calib, axis=1)[:, -n]
    nth_prob_mc_calib = np.sort(mc_calib, axis=1)[:, -n]

    data_counts_uncalib, _ = np.histogram(nth_prob_data_uncalib, bins=bins) 
    mc_counts_uncalib, _   = np.histogram(nth_prob_mc_uncalib, bins=bins)
    data_counts_calib, _ = np.histogram(nth_prob_data_calib, bins=bins) 
    mc_counts_calib, _   = np.histogram(nth_prob_mc_calib, bins=bins)

    ratio_uncalib = np.divide(data_counts_uncalib, mc_counts_uncalib, out=np.zeros_like(mc_counts_uncalib, dtype=float), where=data_counts_uncalib!=0)
    ratio_calib = np.divide(data_counts_calib, mc_counts_calib, out=np.zeros_like(mc_counts_calib, dtype=float), where=data_counts_calib!=0)

    ratio_uncalib=ratio_uncalib*np.sum(mc_counts_uncalib)/np.sum(data_counts_uncalib)
    ratio_calib=ratio_calib*np.sum(mc_counts_calib)/np.sum(data_counts_calib)

    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    ax_ratio.errorbar(bin_centers, ratio_uncalib, yerr=0, fmt='-', drawstyle='steps-mid', markersize=4, color='red', label='Uncalibrated Data/MC')
    ax_ratio.errorbar(bin_centers, ratio_calib, yerr=0, fmt='-', drawstyle='steps-mid', markersize=4, color='blue', label='Calibrated Data/MC')
    ax_ratio.legend()
    ax_ratio.axhline(1, color='gray', linestyle='--', linewidth=1)
    if n==1:
        ax_ratio.axvline(0.4, color='gray', linestyle='--', linewidth=1)

    ratio_data = np.divide(data_counts_calib, data_counts_uncalib, out=np.zeros_like(data_counts_uncalib, dtype=float), where=data_counts_uncalib!=0)
    ratio_mc = np.divide(mc_counts_calib, mc_counts_uncalib, out=np.zeros_like(mc_counts_uncalib, dtype=float), where=data_counts_calib!=0)

    ratio_data=ratio_data*np.sum(data_counts_uncalib)/np.sum(data_counts_calib)
    ratio_mc=ratio_mc*np.sum(mc_counts_uncalib)/np.sum(mc_counts_calib)

    ax_ratio2.errorbar(bin_centers, ratio_data, yerr=0, fmt='-', drawstyle='steps-mid', markersize=4, color='orange', label='Calibrated/Uncalibrated Data')
    ax_ratio2.errorbar(bin_centers, ratio_mc, yerr=0, fmt='-', drawstyle='steps-mid', markersize=4, color='purple', label='Calibrated/Uncalibrated MC')
    ax_ratio2.legend()
    ax_ratio2.axhline(1, color='gray', linestyle='--', linewidth=1)
    if n==1:
        ax_ratio2.axvline(0.4, color='gray', linestyle='--', linewidth=1)
    ax_ratio2.set_ylim(0, 2)

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        plt.savefig(save_filename, bbox_inches='tight', dpi=300)

def plot_cdf_data_vs_mc(mc_uncalib, mc_calib, data_uncalib, data_calib, n=1, num_bins=10, save_filename=None, title=None):
    
    cdf_ax = plot_prob_cdf(mc_uncalib, n=n, num_bins=num_bins, color='#d2c271', label="MC Uncalibrated")
    cdf_ax = plot_prob_cdf(mc_calib, n=n, num_bins=num_bins, color='#CA5454', ax = cdf_ax, save_filename=save_filename, label="MC Calibrated")
    cdf_ax = plot_prob_cdf(data_uncalib, n=n, num_bins=num_bins, color='#40B0A6', ax = cdf_ax, label="Data Uncalibrated")
    cdf_ax = plot_prob_cdf(data_calib, n=n, num_bins=num_bins, color='#9285d5', ax = cdf_ax, save_filename=save_filename, label="Data Calibrated")
    cdf_ax.set_xlim(0, 1/n)

    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n if n < 20 else n % 10, 'th')
    if title is None:
        cdf_ax.set_title(f'Probability Distribution of the {n}{suffix} Most Probable Plane')

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        plt.savefig(save_filename, bbox_inches='tight', dpi=300)


def plot_prob_cdf(probs, n=1, num_bins=10, ax=None, label=None, color=None, save_filename=None, title=None):
    """
    Plots a histogram of the plane numbers that correspond to the n-th highest 
    probability for each event and optionally saves it to a file.

    Parameters:
    - probs: (numpy.ndarray) Matrix of shape [N_events, N_classes] containing probabilities.
    - n: (int) The probability to plot. 1 = Most Probable, 2 = 2nd Most Probable, etc.
    - num_classes: (int) Total number of planes. If None, inferred from probs.
    - ax: (matplotlib.axes.Axes) Optional axis to plot on. If None, creates a new figure.
    - label: (str) Label for the legend.
    - color: (str) Color of the histogram line.
    - save_filename: (str) Full path/filename to save the plot (e.g., 'rank1_calibrated.pdf').
                     If None, the plot is not saved automatically.
    """
    
    # Sort and pick the n-th highest probability (-1 is highest, -2 is 2nd highest, etc.)
    nth_prob = np.sort(probs, axis=1)[:, -n]

    # Create figure if axis not provided otherwise use existing figure, this allows drawing multiple plots overlaid, if ever desired
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.linspace(0, 1, num_bins + 1)
    num_evts = len(nth_prob)
    weights = np.full(num_evts, 1/num_evts)
    ax.hist(nth_prob, bins=bins, histtype='step', cumulative=True, linewidth=2, label=label, color=color, weights=weights)
    ax.set_ylim(0, 1)
    ax.set_xlabel('CNN Probability Score')
    ax.set_ylabel(f'Number of Events')
    if title is None:
        ax.set_title('Event Probability Distribtuion')
    else:
        ax.set_title(title)
    
    if label:
        ax.legend()

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        ax.figure.savefig(save_filename, bbox_inches='tight', dpi=300)
        
    return ax

#BEGIN - WORK IN PROGESS
# Below are early steps taken to implement a diagonal Dirichlet calibration as an alternative to SMS calibration

# A local,  version of DiagonalDirichletCalibrator since versioning issues mean we're stuck with an older version with a known bug
# so we're here implementing the fixed version here ourselves
# In future, if we get access to a newer python version, enabling an upgrade to a newer dirichletcal version, this may be removed
class DiagonalDirichletCalibrator(dirichletcal.calib.fulldirichlet.FullDirichletCalibrator):
    def fit(self, X, y, X_val=None, y_val=None, *args, **kwargs):

        self.weights_ = self.weights_init

        if X_val is None:
            X_val = X.copy()
            y_val = y.copy()

        _X = np.copy(X)
        _X = np.log(dirichletcal.utils.clip_for_log(_X))
        _X_val = np.copy(X_val)
        _X_val = np.log(dirichletcal.utils.clip_for_log(X_val))

        self.calibrator_ = dirichletcal.calib.multinomial.MultinomialRegression(
            method='Diag', reg_lambda=self.reg_lambda, reg_mu=self.reg_mu,
            reg_norm=self.reg_norm, ref_row=self.ref_row,
            optimizer=self.optimizer)
        print(f"Probability matrix shape: {_X.shape}")
        print(f"y shape: {y.shape}")
        self.calibrator_.fit(_X, y, labels = np.arange(0, 215), *args, **kwargs)
        self.final_loss = log_loss(y_val, self.calibrator_.predict_proba(_X_val))

        return self

# Manual extension of the probmetrics or probmetrics style classes for diagonal Dirichlet calibration
class DiagDirichletCalibrator(Calibrator):
    def __init__(
        self,
        n_cv: int = 0,
        use_odir: bool = False,
        reg_lambda: float = 0.0,
        reg_mu: Optional[float] = None,
        reg_lambda_grid: Optional[List[float]] = None,
        reg_mu_grid: Optional[List[float]] = None,
    ):
        super().__init__()

        self.n_cv = n_cv
        self.use_odir = use_odir
        self.reg_lambda = reg_lambda
        self.reg_mu = reg_mu
        self.reg_lambda_grid = reg_lambda_grid
        self.reg_mu_grid = reg_mu_grid

    def _fit_impl(self, X: np.ndarray, y: np.ndarray) -> None:
        self.cal_ = DiagonalDirichletCalibrator()
        self.cal_.l2 = self.reg_lambda
        self.cal_.fit(X, y)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        return self.cal_.predict_proba(X)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="SMS Calibrator. This scirpt is dual purpose. When no model is supplied, it calibrates it's own. And it also produces a number of relevant plots")
    parser.add_argument("-i", "--input", required=True, type=str, help="Path to text file containing raw DANN output")
    parser.add_argument("-m", "--model", required=False, type=str, help="Path to pre-existing model")
    parser.add_argument("-o", "--output", default="reliability_diagram.png", type=str, help="Output root file")
    parser.add_argument("-d", "--data", type=str, help="Input data file")
    parser.add_argument("-s", "--samplerate", type=float, help="Proportion of input files to sample")
    parser.add_argument("--split_data",action='store_true', help="Flag determining if we split our data into training and testing")
    parser.add_argument("-c", "--calibrator", type=str, help="Calibrator. Recommended options: sms or dircal")

    args = parser.parse_args()
    file_path = args.input
    model = args.model
    datafile = args.data
    samplerate = 0.1 
    calibratorstr = args.calibrator
    if args.samplerate is not None:
        samplerate = args.samplerate
    
    justTrain = False #flag to determine if we are only training and not testing on the supplied data
    #Try to load data
    try:
        X_probs, y = load_data_buffered(file_path, samplerate)
    except Exception as e:
        print(f"Error loading data: {e}")
        exit()

    #print(f"Label range: {y.min()} to {y.max()}")
    print(f"Probability matrix shape: {X_probs.shape}")

    # 5. Standard SMS Calibration
    print("\n--- Starting Calibration ---")
    print("Initializing Probmetrics Calibrator...")
    

    #If an existing calibrator is provided, use that and don't waste a lot of time retraining a new one
    #Otherwise, train one
    if model is None:
        print("No pretrained calibration model provided, making one using the supplied data")
        #Split Data. Note to self - maybe we shouldn't do this. Maybe we should train on a whole playlist and then apply it to other playlists
        if args.split_data:
            print("Splitting data into Train and Test sets...")
            X_train_probs, X_test_probs, y_train, y_test = train_test_split(
            X_probs, y, test_size=0.2, random_state=42
            )
            del X_probs, y # Free up memory (delete the full original arrays which can be quite big and eat a lot of RAM)
        else:
            print("You haven't provided a premade calibrator and also not set the --split-data flag to split your data into test/train sets. So I will assume that you don't wish to test on this data, and instead train on all of the supplied data.")
            justTrain = True
            X_train_probs = X_probs
            y_train = y
        #print("X_test_probs: ", X_test_probs)
        #print("y_test: ", y_test)
        #print("y_train: ", y_train)
        #print(f"X test max: {max(X_test_probs[0])}")
        
        if calibratorstr is None:
            calibrator = SMSCalibrator (
            penalty = 'ridge' ,
            rho = 1.0 ,
            tau = 1.0 ,
            lambda_intercept = 1.0 ,
            lambda_diagonal = 1.0 ,
            lambda_off_diagonal = 1.0)

            calibratorstr = "sms"
        elif calibratorstr == "dircal":
            calibrator = DiagDirichletCalibrator(
                reg_lambda = 1.0
            )
            """calibrator = DiagDirichletCalibrator (
            n_cv = 0,
            use_odir = True,
            reg_lambda = 0.0,
            reg_mu = None,
            reg_lambda_grid = None,
            reg_mu_grid = None) """
            print("Setting ODIR for diriclet calibrator")
        else:
            calibrator = get_calibrator(calibratorstr)


        print("==============================================")
        print(f"Using {calibratorstr} calibrator")
        print("==============================================")
        print(f"Fitting Calibrator on {len(y_train)} samples (this may take time)...")
        print(f"X_train_probs: {X_train_probs.shape}")
        print(f"y_train: {y_train.shape}")
        calibrator.fit(X_train_probs, y_train)
        if args.split_data:
            test_probs = X_test_probs
            y_probs = y_test
    else:
        print(f"Loading model from {model}...")
        calibrator = joblib.load(model)
        test_probs = X_probs
        y_probs = y
        calibratorstr = model.split("/")[-1]
        calibratorstr = calibratorstr.split("_")[0]
        print("==============================================")
        print(f"Using premade {calibratorstr} calibrator")
        print("==============================================")


    #Testing and plotting
    if not justTrain:
        print("Calibrating Test Data...")
        import matplotlib.pyplot as plt #We only need plotting functionality in this code block so lazy-loading it here rather than at the top
        test_calibrated = calibrator.predict_proba(test_probs)

        #Plot uncalibrated reliability
        plot_multirank_reliability(
            test_probs, 
            y_probs, 
            title="Uncalibrated DANN Output", 
            filename=f"{calibratorstr}_reliability_uncalibrated.pdf",
            n_ranks=5
        )

        #Plot calibrated reliability
        plot_multirank_reliability(
            test_calibrated, 
            y_probs, 
            title="SMS Calibrated Output", 
            filename=f"{calibratorstr}_reliability_calibrated.pdf",
            n_ranks=5
        )
        #Plot uncalibrated vs calibrated plane distributions
        n1_ax = plot_nth_plane_distribution(test_probs, color="red", label="Uncalibrated", save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_top_plane.pdf"), title=f'Distribution of the most probable Plane Before Calibration')
        n2_ax = plot_nth_plane_distribution(test_probs, color="red", label="Uncalibrated", n=2, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_second_plane.pdf"), title=f'Distribution of the second most probable Plane Before Calibration')
        n3_ax = plot_nth_plane_distribution(test_probs, color="red", label="Uncalibrated", n=3, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_third_plane.pdf"), title=f'Distribution of the third most probable Plane Before Calibration')

        #Plot calibrated plane distributions

        plot_nth_plane_distribution(test_calibrated, color="blue", label="Calibrated", ax=n1_ax, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_top_plane.pdf"), title=f'Distribution of the most probable Plane')
        plot_nth_plane_distribution(test_calibrated, color="blue", label="Calibrated", ax=n2_ax, n=2, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_second_plane.pdf"), title=f'Distribution of the second most probable Plane')
        plot_nth_plane_distribution(test_calibrated, color="blue", label="Calibrated", ax=n3_ax, n=3, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_third_plane.pdf"), title=f'Distribution of the third most probable Plane')

        #Plot calibrated and uncalibrated probability distributions for most probable plane (not POT normalised)
        cdf_ax = plot_prob_cdf(test_probs, n=1, num_bins=100, color='#d2c271', label="MC Uncalibrated")
        cdf_ax = plot_prob_cdf(test_calibrated, n=1, num_bins=100, color='#CA5454', ax = cdf_ax, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_cdf_n1.pdf"), label="MC Calibrated")

        if datafile is None:
            ax = plot_prob_distributions(test_probs, n=1, num_bins=100, color='#d2c271', save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_prob_distribution.pdf"), label="MC Uncalibrated")
            ax = plot_prob_distributions(test_calibrated, n=1, num_bins=100, ax=ax, color='#CA5454', save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_prob_distribution.pdf"), label="MC Calibrated")
        else:
            X_probs_data, y_none = load_data_buffered(datafile, samplerate)
            test_calibrated_data = calibrator.predict_proba(X_probs_data)
            plot_prob_distributions_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=1, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_prob_distribution_n1.pdf"))
            plot_prob_distributions_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=2, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_prob_distribution_n2.pdf"))
            plot_prob_distributions_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=3, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_prob_distribution_n3.pdf"))
            plot_cdf_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=1, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_cdf_n1.pdf"))
            plot_cdf_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=2, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_cdf_n2.pdf"))
            plot_cdf_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=3, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_cdf_n3.pdf"))
    print(f"\nProcessing Complete. Results saved to: {os.path.abspath(OUTPUT_DIR)}")
    
    if model is None:
        SAVE_PATH = os.path.join(OUTPUT_DIR, f"{calibratorstr}_calibrator.pkl")
        print(f"Saving calibrated model to {SAVE_PATH}...")

        # This saves the entire calibrator object, including the learned matrix W
        joblib.dump(calibrator, SAVE_PATH)

        print("Model saved successfully.")