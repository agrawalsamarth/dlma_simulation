import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error

# --- 1. PRE-PROCESS DATA (Outside the Loop) ---
# Assuming 'df' is your options data and 'spot_history' is your full price series
# spot_history should have columns ['date', 'price']

# Sort by date to allow binary search
spot_history = spot_history.sort_values('date').reset_index(drop=True)
full_price_values = spot_history['price'].values
full_date_values = spot_history['date'].values

# Calculate integer indices for every option once
all_metadata = []
for row in df.itertuples():
    start_idx = np.searchsorted(full_date_values, row.date_start)
    end_idx = np.searchsorted(full_date_values, row.date_end, side='right')
    
    all_metadata.append({
        'start_idx': start_idx,
        'end_idx': end_idx,
        'strike': row.strike,
        'date_to_exp': row.date_to_exp
    })

y_observed = df['market_price'].values

# --- 2. THE CALIBRATION WRAPPER ---
def calibration_wrapper(combined_input, mu, sigma):
    """
    Inputs:
        combined_input: Tuple of (price_array, metadata_list)
        mu, sigma: Parameters being tuned by curve_fit
    """
    price_array, metadata = combined_input
    predictions = []
    
    for row in metadata:
        # Get a VIEW of the prices (Zero memory overhead)
        relevant_prices = price_array[row['start_idx'] : row['end_idx']]
        
        if len(relevant_prices) == 0:
            predictions.append(0.0)
            continue
            
        # --- Normalization Logic ---
        # Adjust this based on your specific normalization needs
        norm_factor = relevant_prices[0] 
        if norm_factor == 0: norm_factor = 1e-8 # Prevent division by zero
        
        normalized_path = relevant_prices / norm_factor
        
        # Call your fixed pricer
        # Note: If strike needs normalization too, use row['strike'] / norm_factor
        p = option_pricer(
            prices=normalized_path,
            strike=row['strike'],
            date_to_exp=row['date_to_exp'],
            mu=mu,
            sigma=sigma
        )
        predictions.append(p)
        
    return np.array(predictions)

# --- 3. K-FOLD CROSS VALIDATION ---
kf = KFold(n_splits=5, shuffle=True, random_state=42)
results_log = []

print(f"{'Fold':<6} | {'Mu':<10} | {'Sigma':<10} | {'Test RMSE':<10}")
print("-" * 45)

for fold, (train_idx, test_idx) in enumerate(kf.split(y_observed), 1):
    # Split Metadata
    train_meta = [all_metadata[i] for i in train_idx]
    test_meta = [all_metadata[i] for i in test_idx]
    
    # Split Targets
    y_train, y_test = y_observed[train_idx], y_observed[test_idx]
    
    # Pack inputs (Passing reference to full_price_values, not a copy)
    train_inputs = (full_price_values, train_meta)
    test_inputs = (full_price_values, test_meta)
    
    # Fit (In-Sample)
    # Bounds: mu can be anything, sigma must be > 0
    try:
        popt, _ = curve_fit(
            calibration_wrapper, 
            train_inputs, 
            y_train, 
            p0=[0.05, 0.2],
            bounds=([-np.inf, 1e-6], [np.inf, np.inf])
        )
        
        mu_fit, sigma_fit = popt
        
        # Predict (Out-of-Sample)
        y_pred = calibration_wrapper(test_inputs, *popt)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))
        
        results_log.append({'mu': mu_fit, 'sigma': sigma_fit, 'rmse': rmse})
        print(f"{fold:<6} | {mu_fit:<10.4f} | {sigma_fit:<10.4f} | {rmse:<10.4f}")
        
    except Exception as e:
        print(f"Fold {fold} failed: {e}")

# --- 4. FINAL SUMMARY ---
avg_rmse = np.mean([r['rmse'] for r in results_log])
print("-" * 45)
print(f"Average Out-of-Sample RMSE: {avg_rmse:.4f}")
