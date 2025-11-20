# Soil simulation support functions

def infill_soil_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    More robust version that handles missing soil texture data at the horizon level
    rather than excluding entire groups. Works on entire DataFrame.
    
    Args:
        df: Input soil DataFrame
        
    Returns:
        DataFrame with infilled texture data
    """
    # Step 1: Group by 'compname_grp' and process each group
    grouped = df.groupby("compname_grp")
    
    def process_group(group: pd.DataFrame) -> pd.DataFrame:
        """Process each group individually with horizon-level filtering"""
        group = group.copy()
        
        # Identify horizons within top 50cm that have complete texture data
        top_50_mask = group["hzdepb_r"] <= 50
        texture_cols = ["sandtotal_r", "claytotal_r", "silttotal_r"]
        
        # Check which horizons have missing texture data
        missing_texture = group[texture_cols].isnull().any(axis=1)
        
        # For horizons in top 50cm with missing texture data, try to infill or mark for exclusion
        problematic_horizons = top_50_mask & missing_texture
        
        if problematic_horizons.any():
            # Try to infill missing texture data using group averages or nearby horizons
            group = infill_missing_texture_data(group, problematic_horizons)
            
            # If infilling fails, mark those specific horizons as problematic
            # but keep the rest of the group
            still_missing = group[texture_cols].isnull().any(axis=1)
            if (top_50_mask & still_missing).any():
                # Add a flag column to identify problematic horizons
                group['texture_data_complete'] = ~(top_50_mask & still_missing)
            else:
                group['texture_data_complete'] = True
        else:
            group['texture_data_complete'] = True
        
        return group
    
    # Apply processing to all groups and combine results
    processed_groups = []
    for name, group in grouped:
        processed_group = process_group(group)
        processed_groups.append(processed_group)
    
    # Combine all processed groups
    result_df = pd.concat(processed_groups, ignore_index=True)
    
    # Infill the _l and _h values for all valid horizons
    result_df = infill_range_values(result_df)
    
    # Apply RFV imputation if needed
    try:
        result_df = result_df.apply(impute_rfv_values, axis=1)
    except Exception as e:
        print(f"Warning: RFV imputation failed: {e}")
    
    return result_df


# Range Value Infilling Functions

def infill_range_values(df: pd.DataFrame) -> pd.DataFrame:
    """
    Robust infilling of _l and _h values using data-driven approaches
    
    Args:
        df: Input DataFrame
        
    Returns:
        DataFrame with infilled range values
    """
    df = df.copy()
    
    # Define parameter groups with their characteristics
    param_configs = {
        'texture': {
            'params': ['sandtotal', 'claytotal', 'silttotal'],
            'units': '%',
            'sum_constraint': 100,  # Sand + Clay + Silt should sum to ~100%
            'fallback_range': 8
        },
        'bulk_density': {
            'params': ['dbovendry'],
            'units': 'g/cm³',
            'typical_range': (0.8, 2.2),
            'fallback_range': 0.01
        },
        'water_retention': {
            'params': ['wthirdbar', 'wfifteenbar'],
            'units': '%',
            'typical_range': (0, 60),
            'fallback_range': {'wthirdbar': 1, 'wfifteenbar': 0.6}
        }
    }
    
    # Process each parameter group
    for group_name, config in param_configs.items():
        for param in config['params']:
            df = infill_parameter_range(df, param, config, group_name)
    
    # Final validation and constraint enforcement
    df = enforce_range_constraints(df)
    
    return df

def infill_parameter_range(df: pd.DataFrame, param: str, config: dict, group_name: str) -> pd.DataFrame:
    """
    Infill _l and _h values for a specific parameter using multiple strategies
    
    Args:
        df: Input DataFrame
        param: Parameter name
        config: Parameter configuration
        group_name: Parameter group name
        
    Returns:
        DataFrame with infilled range values
    """
    df = df.copy()
    r_col = f"{param}_r"
    l_col = f"{param}_l"
    h_col = f"{param}_h"
    
    if r_col not in df.columns:
        return df
    
    # Strategy 1: Learn from existing complete ranges in the dataset
    learned_ranges = learn_parameter_ranges(df, param, config)
    
    # Strategy 2: Apply horizon-specific and component-group-specific ranges
    context_ranges = get_contextual_ranges(df, param, config)
    
    # Strategy 3: Infill missing _l values
    if l_col not in df.columns:
        df[l_col] = np.nan
    
    missing_l_mask = df[l_col].isnull() & df[r_col].notna()
    if missing_l_mask.any():
        df.loc[missing_l_mask, l_col] = df.loc[missing_l_mask].apply(
            lambda row: calculate_lower_bound(row, param, learned_ranges, context_ranges, config),
            axis=1
        )
    
    # Strategy 4: Infill missing _h values  
    if h_col not in df.columns:
        df[h_col] = np.nan
        
    missing_h_mask = df[h_col].isnull() & df[r_col].notna()
    if missing_h_mask.any():
        df.loc[missing_h_mask, h_col] = df.loc[missing_h_mask].apply(
            lambda row: calculate_upper_bound(row, param, learned_ranges, context_ranges, config),
            axis=1
        )
    
    # Apply bounds checking
    df[l_col] = df[l_col].apply(lambda x: max(x, 0) if pd.notna(x) else x)
    if 'typical_range' in config:
        max_val = config['typical_range'][1]
        df[h_col] = df[h_col].apply(lambda x: min(x, max_val) if pd.notna(x) else x)
    
    return df

def learn_parameter_ranges(df: pd.DataFrame, param: str, config: dict) -> dict:
    """
    Learn typical ranges from existing complete data
    
    Args:
        df: Input DataFrame
        param: Parameter name
        config: Parameter configuration
        
    Returns:
        Dictionary of learned ranges by context
    """
    r_col = f"{param}_r"
    l_col = f"{param}_l"
    h_col = f"{param}_h"
    
    # Get complete records (have all three values)
    complete_mask = df[[r_col, l_col, h_col]].notna().all(axis=1)
    complete_data = df[complete_mask]
    
    if complete_data.empty:
        return {'default_spread': config.get('fallback_range', 5)}
    
    # Calculate actual spreads
    lower_spreads = complete_data[r_col] - complete_data[l_col]
    upper_spreads = complete_data[h_col] - complete_data[r_col]
    
    ranges_by_context = {}
    
    # Learn ranges by horizon type
    if 'hzname' in complete_data.columns:
        for hzname in complete_data['hzname'].dropna().unique():
            hz_data = complete_data[complete_data['hzname'] == hzname]
            if len(hz_data) >= 3:  # Need minimum sample size
                hz_lower = hz_data[r_col] - hz_data[l_col]
                hz_upper = hz_data[h_col] - hz_data[r_col]
                
                ranges_by_context[f'hzname_{hzname}'] = {
                    'lower_spread_median': hz_lower.median(),
                    'upper_spread_median': hz_upper.median(),
                    'lower_spread_std': hz_lower.std(),
                    'upper_spread_std': hz_upper.std(),
                    'sample_size': len(hz_data)
                }
    
    # Learn ranges by component group
    if 'compname_grp' in complete_data.columns:
        for grp in complete_data['compname_grp'].dropna().unique():
            grp_data = complete_data[complete_data['compname_grp'] == grp]
            if len(grp_data) >= 5:
                grp_lower = grp_data[r_col] - grp_data[l_col]
                grp_upper = grp_data[h_col] - grp_data[r_col]
                
                ranges_by_context[f'compgrp_{grp}'] = {
                    'lower_spread_median': grp_lower.median(),
                    'upper_spread_median': grp_upper.median(),
                    'lower_spread_std': grp_lower.std(),
                    'upper_spread_std': grp_upper.std(),
                    'sample_size': len(grp_data)
                }
    
    # Overall dataset statistics
    ranges_by_context['overall'] = {
        'lower_spread_median': lower_spreads.median(),
        'upper_spread_median': upper_spreads.median(),
        'lower_spread_std': lower_spreads.std(),
        'upper_spread_std': upper_spreads.std(),
        'sample_size': len(complete_data)
    }
    
    return ranges_by_context

def get_contextual_ranges(df: pd.DataFrame, param: str, config: dict) -> dict:
    """
    Get parameter-specific contextual ranges based on soil science knowledge
    
    Args:
        df: Input DataFrame
        param: Parameter name
        config: Parameter configuration
        
    Returns:
        Dictionary of contextual ranges
    """
    contextual_ranges = {}
    
    # Texture-specific ranges
    if param in ['sandtotal', 'claytotal', 'silttotal']:
        contextual_ranges.update({
            'sandtotal': {'typical_spread': 12, 'min_spread': 5, 'max_spread': 25},
            'claytotal': {'typical_spread': 8, 'min_spread': 3, 'max_spread': 20},
            'silttotal': {'typical_spread': 10, 'min_spread': 4, 'max_spread': 22}
        })
    
    # Bulk density ranges by horizon type
    elif param == 'dbovendry':
        contextual_ranges.update({
            'A_horizons': {'typical_spread': 0.15, 'range': (0.8, 1.6)},
            'B_horizons': {'typical_spread': 0.12, 'range': (1.0, 1.8)},
            'C_horizons': {'typical_spread': 0.10, 'range': (1.2, 2.0)}
        })
    
    # Water retention ranges
    elif param in ['wthirdbar', 'wfifteenbar']:
        contextual_ranges.update({
            'wthirdbar': {'typical_spread': 3, 'clay_factor': 0.4},
            'wfifteenbar': {'typical_spread': 2, 'clay_factor': 0.3}
        })
    
    return contextual_ranges

def calculate_lower_bound(row: pd.Series, param: str, learned_ranges: dict, 
                         context_ranges: dict, config: dict) -> float:
    """
    Calculate appropriate lower bound for a parameter
    
    Args:
        row: DataFrame row
        param: Parameter name
        learned_ranges: Learned ranges from data
        context_ranges: Contextual ranges
        config: Parameter configuration
        
    Returns:
        Calculated lower bound
    """
    r_value = row[f"{param}_r"]
    if pd.isna(r_value):
        return np.nan
    
    # Priority order for range estimation
    spread_estimates = []
    
    # 1. Try horizon-specific learned range
    if 'hzname' in row and pd.notna(row['hzname']):
        hz_key = f"hzname_{row['hzname']}"
        if hz_key in learned_ranges and learned_ranges[hz_key]['sample_size'] >= 3:
            spread = learned_ranges[hz_key]['lower_spread_median']
            if pd.notna(spread) and spread > 0:
                spread_estimates.append(('horizon_learned', spread))
    
    # 2. Try component group learned range
    if 'compname_grp' in row and pd.notna(row['compname_grp']):
        grp_key = f"compgrp_{row['compname_grp']}"
        if grp_key in learned_ranges and learned_ranges[grp_key]['sample_size'] >= 5:
            spread = learned_ranges[grp_key]['lower_spread_median']
            if pd.notna(spread) and spread > 0:
                spread_estimates.append(('group_learned', spread))
    
    # 3. Use overall dataset learned range
    if 'overall' in learned_ranges:
        spread = learned_ranges['overall']['lower_spread_median']
        if pd.notna(spread) and spread > 0:
            spread_estimates.append(('overall_learned', spread))
    
    # 4. Use contextual/pedological knowledge
    if param in context_ranges:
        ctx_spread = context_ranges[param].get('typical_spread', config.get('fallback_range', 5))
        
        # Adjust spread based on soil properties
        if param in ['wthirdbar', 'wfifteenbar'] and 'claytotal_r' in row:
            clay_content = row['claytotal_r']
            if pd.notna(clay_content):
                clay_factor = context_ranges[param].get('clay_factor', 0.3)
                ctx_spread = ctx_spread + (clay_content * clay_factor / 100)
        
        spread_estimates.append(('contextual', ctx_spread))
    
    # 5. Fallback to default
    if not spread_estimates:
        fallback = config.get('fallback_range', 5)
        spread_estimates.append(('fallback', fallback))
    
    # Use the first (highest priority) estimate
    spread_method, spread = spread_estimates[0]
    
    # Calculate lower bound
    lower_bound = r_value - spread
    
    # Apply parameter-specific constraints
    lower_bound = max(lower_bound, 0)  # Never below 0
    
    return lower_bound

def calculate_upper_bound(row: pd.Series, param: str, learned_ranges: dict, 
                         context_ranges: dict, config: dict) -> float:
    """
    Calculate appropriate upper bound for a parameter
    
    Args:
        row: DataFrame row
        param: Parameter name
        learned_ranges: Learned ranges from data
        context_ranges: Contextual ranges
        config: Parameter configuration
        
    Returns:
        Calculated upper bound
    """
    r_value = row[f"{param}_r"]
    if pd.isna(r_value):
        return np.nan
    
    # Similar logic to lower bound but for upper spreads
    spread_estimates = []
    
    # 1. Horizon-specific learned range
    if 'hzname' in row and pd.notna(row['hzname']):
        hz_key = f"hzname_{row['hzname']}"
        if hz_key in learned_ranges and learned_ranges[hz_key]['sample_size'] >= 3:
            spread = learned_ranges[hz_key]['upper_spread_median']
            if pd.notna(spread) and spread > 0:
                spread_estimates.append(('horizon_learned', spread))
    
    # 2. Component group learned range
    if 'compname_grp' in row and pd.notna(row['compname_grp']):
        grp_key = f"compgrp_{row['compname_grp']}"
        if grp_key in learned_ranges and learned_ranges[grp_key]['sample_size'] >= 5:
            spread = learned_ranges[grp_key]['upper_spread_median']
            if pd.notna(spread) and spread > 0:
                spread_estimates.append(('group_learned', spread))
    
    # 3. Overall dataset learned range
    if 'overall' in learned_ranges:
        spread = learned_ranges['overall']['upper_spread_median']
        if pd.notna(spread) and spread > 0:
            spread_estimates.append(('overall_learned', spread))
    
    # 4. Contextual knowledge
    if param in context_ranges:
        ctx_spread = context_ranges[param].get('typical_spread', config.get('fallback_range', 5))
        
        # Adjust for soil properties
        if param in ['wthirdbar', 'wfifteenbar'] and 'claytotal_r' in row:
            clay_content = row['claytotal_r']
            if pd.notna(clay_content):
                clay_factor = context_ranges[param].get('clay_factor', 0.3)
                ctx_spread = ctx_spread + (clay_content * clay_factor / 100)
        
        spread_estimates.append(('contextual', ctx_spread))
    
    # 5. Fallback
    if not spread_estimates:
        fallback = config.get('fallback_range', 5)
        spread_estimates.append(('fallback', fallback))
    
    spread_method, spread = spread_estimates[0]
    upper_bound = r_value + spread
    
    # Apply parameter-specific upper limits
    if 'typical_range' in config:
        upper_bound = min(upper_bound, config['typical_range'][1])
    
    return upper_bound

def enforce_range_constraints(df: pd.DataFrame) -> pd.DataFrame:
    """
    Enforce logical constraints on range values
    
    Args:
        df: Input DataFrame
        
    Returns:
        DataFrame with enforced constraints
    """
    df = df.copy()
    
    # Texture constraints: ensure _l ≤ _r ≤ _h and sum constraints
    texture_params = ['sandtotal', 'claytotal', 'silttotal']
    
    for param in texture_params:
        r_col, l_col, h_col = f"{param}_r", f"{param}_l", f"{param}_h"
        
        if all(col in df.columns for col in [r_col, l_col, h_col]):
            # Ensure _l ≤ _r ≤ _h
            df[l_col] = np.minimum(df[l_col], df[r_col])
            df[h_col] = np.maximum(df[h_col], df[r_col])
    
    # Check texture sum constraints (sand + clay + silt ≈ 100%)
    has_all_texture = all(f"{p}_r" in df.columns for p in texture_params)
    if has_all_texture:
        texture_sum = df[['sandtotal_r', 'claytotal_r', 'silttotal_r']].sum(axis=1)
        # Flag records where sum is far from 100%
        sum_issues = ((texture_sum < 95) | (texture_sum > 105)) & texture_sum.notna()
        if sum_issues.any():
            print(f"Warning: {sum_issues.sum()} records have texture sums outside 95-105% range")
    
    # Water retention constraints: wfifteenbar ≤ wthirdbar
    if 'wfifteenbar_r' in df.columns and 'wthirdbar_r' in df.columns:
        # Fifteen bar should be ≤ third bar (higher pressure = less water)
        invalid_water = df['wfifteenbar_r'] > df['wthirdbar_r']
        if invalid_water.any():
            print(f"Warning: {invalid_water.sum()} records have wfifteenbar > wthirdbar")
    
    return df

# Texture Data Recovery Functions
def infill_missing_texture_data(group: pd.DataFrame, problematic_mask: pd.Series) -> pd.DataFrame:
    """
    Attempt to infill missing texture data using pedologically-informed strategies
    
    Args:
        group: DataFrame group to process
        problematic_mask: Boolean mask for problematic horizons
        
    Returns:
        DataFrame with infilled texture data
    """
    group = group.copy()
    texture_cols = ["sandtotal_r", "claytotal_r", "silttotal_r"]
    
    # Ensure we have depth information
    if 'hzdept_r' not in group.columns:
        print("Warning: Missing hzdept_r (top depth), using 0 for surface horizons")
        group['hzdept_r'] = group.get('hzdept_r', 0)
    
    # Strategy 1: PRIORITY - Match by horizon name (pedologically similar horizons)
    if 'hzname' in group.columns:
        group = horizon_name_texture_infill(group, texture_cols, problematic_mask)
    
    # Check if we still have missing data after horizon name matching
    still_missing = group[texture_cols].isnull().any(axis=1)
    if not still_missing.any():
        return group  # All data recovered, no need for further strategies
    
    # Strategy 2: Depth-weighted averaging from similar depth ranges
    group = depth_weighted_infill(group, texture_cols, problematic_mask)
    
    # Check progress after depth-weighted infill
    still_missing = group[texture_cols].isnull().any(axis=1)
    if not still_missing.any():
        return group
    
    # Strategy 3: Vertical interpolation within individual soil components
    group = within_component_interpolation(group, texture_cols)
    
    # Check progress after interpolation
    still_missing = group[texture_cols].isnull().any(axis=1)
    if not still_missing.any():
        return group
    
    # Strategy 4: Cross-component interpolation for similar depths
    group = cross_component_depth_interpolation(group, texture_cols)
    
    # Strategy 5: Fallback to group statistics for any remaining missing values
    for col in texture_cols:
        if group[col].isnull().any():
            # Use depth-weighted group mean as last resort
            depth_weighted_mean = calculate_depth_weighted_mean(group, col)
            if not pd.isna(depth_weighted_mean):
                group[col] = group[col].fillna(depth_weighted_mean)
    
    return group

def horizon_name_texture_infill(group: pd.DataFrame, texture_cols: List[str], problematic_mask: pd.Series) -> pd.DataFrame:
    """
    Infill missing texture data using horizons with matching or similar names
    
    Args:
        group: DataFrame group to process
        texture_cols: List of texture column names
        problematic_mask: Boolean mask for problematic horizons
        
    Returns:
        DataFrame with infilled values
    """
    group = group.copy()
    
    for idx in group[problematic_mask].index:
        if idx not in group.index:
            continue
        
        target_hzname = group.loc[idx, 'hzname']
        if pd.isna(target_hzname):
            continue
            
        # Clean and standardize horizon name for matching
        target_hzname_clean = standardize_horizon_name(target_hzname)
        
        for col in texture_cols:
            if pd.isna(group.loc[idx, col]):
                # Find matching horizon names with valid data
                matching_values = []
                similarity_scores = []
                
                for other_idx in group.index:
                    if other_idx == idx or pd.isna(group.loc[other_idx, col]):
                        continue
                    
                    other_hzname = group.loc[other_idx, 'hzname']
                    if pd.isna(other_hzname):
                        continue
                    
                    other_hzname_clean = standardize_horizon_name(other_hzname)
                    similarity = calculate_horizon_similarity(target_hzname_clean, other_hzname_clean)
                    
                    if similarity > 0.5:  # Only use reasonably similar horizons
                        matching_values.append(group.loc[other_idx, col])
                        similarity_scores.append(similarity)
                
                if matching_values and similarity_scores:
                    # Use similarity-weighted average
                    if len(matching_values) == 1:
                        infilled_value = matching_values[0]
                    else:
                        infilled_value = np.average(matching_values, weights=similarity_scores)
                    
                    group.loc[idx, col] = infilled_value
                    
                    # Track infilling method for reporting
                    if 'infill_method' not in group.columns:
                        group['infill_method'] = ''
                    group.loc[idx, 'infill_method'] += f'{col}:hzname({target_hzname_clean}); '
    
    return group

def standardize_horizon_name(hzname: str) -> str:
    """
    Standardize horizon names for better matching
    
    Args:
        hzname: Raw horizon name
        
    Returns:
        Cleaned and standardized horizon name
    """
    if pd.isna(hzname):
        return ''
    
    # Convert to string and clean
    hzname = str(hzname).strip().upper()
    
    # Remove common suffixes/prefixes that don't affect texture matching
    # Remove numbers at the end (layer designations)
    hzname = re.sub(r'\d+$', '', hzname)
    
    # Remove common punctuation but keep important ones like '/'
    hzname = re.sub(r'[^\w/]', '', hzname)
    
    return hzname

def calculate_horizon_similarity(hz1: str, hz2: str) -> float:
    """
    Calculate similarity between two horizon names
    
    Args:
        hz1, hz2: Horizon names to compare
        
    Returns:
        Similarity score from 0 (no similarity) to 1 (identical)
    """
    if hz1 == hz2:
        return 1.0
    
    if not hz1 or not hz2:
        return 0.0
    
    # Check for main horizon letter match
    main_hz1 = hz1[0] if hz1 else ''
    main_hz2 = hz2[0] if hz2 else ''
    
    if main_hz1 == main_hz2:
        base_score = 0.8  # Same main horizon type
        
        # Bonus for additional character matches
        common_chars = set(hz1) & set(hz2)
        bonus = len(common_chars) / max(len(hz1), len(hz2)) * 0.2
        
        return min(1.0, base_score + bonus)
    
    # Check for related horizons (pedologically similar)
    related_groups = [
        {'A', 'AP', 'AE'},  # Surface horizons
        {'E', 'EB', 'BE'},  # Eluvial horizons
        {'B', 'BT', 'BW', 'BC', 'BS'},  # Illuvial/subsurface horizons
        {'C', 'CB', 'CR'},  # Parent material
        {'O', 'OA', 'OE'},  # Organic horizons
    ]
    
    for group in related_groups:
        if main_hz1 in group and main_hz2 in group:
            return 0.6  # Related horizon types
    
    # Handle transitional horizons like 'B/C', 'A/E'
    if '/' in hz1 or '/' in hz2:
        parts1 = hz1.split('/') if '/' in hz1 else [hz1]
        parts2 = hz2.split('/') if '/' in hz2 else [hz2]
        
        max_sim = 0
        for p1 in parts1:
            for p2 in parts2:
                sim = calculate_horizon_similarity(p1.strip(), p2.strip())
                max_sim = max(max_sim, sim)
        
        return max_sim * 0.8  # Slight penalty for composite matching
    
    return 0.0  # No similarity

def depth_weighted_infill(group: pd.DataFrame, texture_cols: List[str], problematic_mask: pd.Series) -> pd.DataFrame:
    """
    Fill missing values using depth-weighted averages from similar depth ranges
    
    Args:
        group: DataFrame group to process
        texture_cols: List of texture column names
        problematic_mask: Boolean mask for problematic horizons
        
    Returns:
        DataFrame with infilled values
    """
    group = group.copy()
    
    for idx in group[problematic_mask].index:
        if idx not in group.index:
            continue
            
        target_top = group.loc[idx, 'hzdept_r']
        target_bottom = group.loc[idx, 'hzdepb_r']
        target_mid = (target_top + target_bottom) / 2
        
        # Find horizons with overlapping or similar depth ranges
        for col in texture_cols:
            if pd.isna(group.loc[idx, col]):
                similar_depth_values = []
                weights = []
                
                for other_idx in group.index:
                    if other_idx == idx or pd.isna(group.loc[other_idx, col]):
                        continue
                    
                    other_top = group.loc[other_idx, 'hzdept_r']
                    other_bottom = group.loc[other_idx, 'hzdepb_r']
                    other_mid = (other_top + other_bottom) / 2
                    
                    # Calculate overlap and depth similarity
                    overlap = max(0, min(target_bottom, other_bottom) - max(target_top, other_top))
                    depth_similarity = 1 / (1 + abs(target_mid - other_mid))
                    
                    if overlap > 0 or abs(target_mid - other_mid) < 20:  # 20cm tolerance
                        weight = overlap + depth_similarity
                        similar_depth_values.append(group.loc[other_idx, col])
                        weights.append(weight)
                
                if similar_depth_values and weights:
                    weighted_avg = np.average(similar_depth_values, weights=weights)
                    group.loc[idx, col] = weighted_avg
    
    return group

def within_component_interpolation(group: pd.DataFrame, texture_cols: List[str]) -> pd.DataFrame:
    """
    Interpolate missing values within individual soil components based on depth
    
    Args:
        group: DataFrame group to process
        texture_cols: List of texture column names
        
    Returns:
        DataFrame with interpolated values
    """
    group = group.copy()
    
    # Group by individual soil component if component ID exists
    if 'cokey' in group.columns:
        component_groups = group.groupby('cokey')
    else:
        # Fallback: treat entire group as one component
        component_groups = [(None, group)]
    
    processed_components = []
    
    for comp_id, comp_data in component_groups:
        comp_data = comp_data.copy().sort_values('hzdepb_r')
        
        for col in texture_cols:
            if comp_data[col].isnull().any() and comp_data[col].notna().any():
                # Interpolate missing values based on depth
                comp_data[col] = comp_data[col].interpolate(method='index')
        
        processed_components.append(comp_data)
    
    if len(processed_components) > 1:
        return pd.concat(processed_components, ignore_index=False)
    else:
        return processed_components[0] if processed_components else group

def cross_component_depth_interpolation(group: pd.DataFrame, texture_cols: List[str]) -> pd.DataFrame:
    """
    Use data from other components at similar depths to fill missing values
    
    Args:
        group: DataFrame group to process
        texture_cols: List of texture column names
        
    Returns:
        DataFrame with infilled values
    """
    group = group.copy()
    
    for col in texture_cols:
        missing_mask = group[col].isnull()
        if not missing_mask.any():
            continue
        
        for idx in group[missing_mask].index:
            target_depth_mid = (group.loc[idx, 'hzdept_r'] + group.loc[idx, 'hzdepb_r']) / 2
            
            # Find values from similar depths in other components
            depth_matches = []
            for other_idx in group.index:
                if other_idx == idx or pd.isna(group.loc[other_idx, col]):
                    continue
                
                other_depth_mid = (group.loc[other_idx, 'hzdept_r'] + group.loc[other_idx, 'hzdepb_r']) / 2
                depth_diff = abs(target_depth_mid - other_depth_mid)
                
                if depth_diff <= 15:  # 15cm tolerance for similar depths
                    weight = 1 / (1 + depth_diff)  # Closer depths get higher weight
                    depth_matches.append((group.loc[other_idx, col], weight))
            
            if depth_matches:
                values, weights = zip(*depth_matches)
                weighted_value = np.average(values, weights=weights)
                group.loc[idx, col] = weighted_value
    
    return group

def calculate_depth_weighted_mean(group: pd.DataFrame, col: str) -> float:
    """
    Calculate depth-weighted mean for a column, giving more weight to thicker horizons
    
    Args:
        group: DataFrame group
        col: Column name
        
    Returns:
        Depth-weighted mean value
    """
    valid_data = group[group[col].notna()]
    if valid_data.empty:
        return np.nan
    
    # Calculate horizon thickness as weight
    thicknesses = valid_data['hzdepb_r'] - valid_data['hzdept_r']
    values = valid_data[col]
    
    if thicknesses.sum() == 0:
        return values.mean()
    
    return np.average(values, weights=thicknesses)

# Function to infill missing soil rock fragment data.
def impute_rfv_values(row):
    """
    Improved RFV (Rock Fragment Volume) imputation with better logic and safety checks
    
    RFV represents the percentage of rock fragments by volume in soil horizons.
    Typical ranges: 0-85% (rarely above 85% in productive soils)
    
    Special case: If rfv_r = 0 OR rfv_r is NA, set to 0.02% with range 0.01-0.03% for simulation compatibility.
    """
    
    # Handle missing RFV or zero RFV - both mean "no significant rock fragments"
    if pd.isna(row.get("rfv_r")) or row.get("rfv_r") == 0:
        row["rfv_r"] = 0.02
        row["rfv_l"] = 0.01
        row["rfv_h"] = 0.03
        return row
    
    rfv_r = row["rfv_r"]
    
    # Ensure rfv_r is within reasonable bounds (0.01-85%)
    rfv_r = max(0.01, min(rfv_r, 85))
    row["rfv_r"] = rfv_r
    
    # Determine appropriate range based on rfv_r value and context
    if 0.01 <= rfv_r <= 5:
        # Low RFV: careful with lower bound to avoid going to 0
        lower_spread = min(rfv_r - 0.01, 2)  # Don't go below 0.01
        upper_spread = 4
        
    elif 5 < rfv_r <= 15:
        # Moderate RFV: standard range
        lower_spread = 3
        upper_spread = 5
        
    elif 15 < rfv_r <= 35:
        # High RFV: wider range reflects more variability
        lower_spread = 5
        upper_spread = 8
        
    elif 35 < rfv_r <= 60:
        # Very high RFV: even wider range
        lower_spread = 8
        upper_spread = 12
        
    else:  # rfv_r > 60
        # Extremely high RFV: largest range but respect upper bound
        lower_spread = 10
        upper_spread = min(15, 85 - rfv_r)  # Don't exceed 85%
    
    # Adjust spreads based on horizon context if available
    if 'hzname' in row and pd.notna(row['hzname']):
        hzname = str(row['hzname']).upper()
        
        # Surface horizons typically more variable due to management
        if hzname.startswith('A') and rfv_r > 0.01:
            upper_spread *= 1.3
            
        # Bt horizons often more uniform
        elif 'BT' in hzname:
            lower_spread *= 0.8
            upper_spread *= 0.8
            
        # C horizons more variable (parent material heterogeneity)
        elif hzname.startswith('C'):
            lower_spread *= 1.2
            upper_spread *= 1.2
    
    # Calculate bounds with safety checks
    rfv_l = max(0.01, rfv_r - lower_spread)  # Never below 0.01% for simulation
    rfv_h = min(85, rfv_r + upper_spread)
    
    # Only update _l and _h if they're missing
    if pd.isna(row.get("rfv_l")):
        row["rfv_l"] = rfv_l
    
    if pd.isna(row.get("rfv_h")):
        row["rfv_h"] = rfv_h
    
    # Final consistency check
    if not pd.isna(row.get("rfv_l")) and not pd.isna(row.get("rfv_h")):
        # Ensure logical order: rfv_l ≤ rfv_r ≤ rfv_h
        row["rfv_l"] = min(row["rfv_l"], rfv_r)
        row["rfv_h"] = max(row["rfv_h"], rfv_r)
        
        # Ensure rfv_l never goes to 0 for simulation compatibility
        if row["rfv_l"] <= 0:
            row["rfv_l"] = 0.01
        
        # Ensure reasonable minimum range for simulation
        if row["rfv_h"] - row["rfv_l"] < 0.02:
            spread = max(0.02, rfv_r * 0.2)  # Minimum 0.02% or 20% of value
            row["rfv_l"] = max(0.01, rfv_r - spread/2)
            row["rfv_h"] = min(85, rfv_r + spread/2)
    
    return row