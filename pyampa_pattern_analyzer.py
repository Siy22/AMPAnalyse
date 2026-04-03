"""
PYAMPA Pattern Analyzer - Structure-Activity Relationship (SAR) Principle Extractor
===================================================================================
This script analyzes PYAMPA mutagenesis data to extract general principles
about amino acid substitutions and their effects on:
- Antimicrobial activity (AMP probability)
- Hemolytic activity (Hemolytic probability)
- Overall fitness

Outputs:
1. Amino acid property analysis
2. Best/worst mutation patterns
3. Position-specific insights
4. General SAR principles
5. Actionable design rules

Author: AMP Design Tool
Date: 2025
For: Final Year Project
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import json

# ============================================================================
# CONFIGURATION
# ============================================================================

PYAMPA_FILE = "YOURSEQ_mutagenesis.csv"
OUTPUT_REPORT = "pyampa_analysis_report.txt"
OUTPUT_PRINCIPLES = "design_principles_YOURSEQ.json"

# Analysis thresholds
GOOD_AMP_THRESHOLD = None    # Will be set to median (data-driven)
GOOD_HEMO_THRESHOLD = None   # Will be set to median (data-driven)
GOOD_FITNESS_THRESHOLD = 0.0 # Fixed at zero (biologically meaningful baseline)

# ============================================================================
# AMINO ACID PROPERTY DATABASE
# ============================================================================

AA_PROPERTIES = {
    # Standard amino acids with their properties
    'A': {'name': 'Alanine', 'type': 'nonpolar', 'size': 'small', 'charge': 'neutral', 'hydrophobic': True, 'aromatic': False},
    'R': {'name': 'Arginine', 'type': 'basic', 'size': 'large', 'charge': 'positive', 'hydrophobic': False, 'aromatic': False},
    'N': {'name': 'Asparagine', 'type': 'polar', 'size': 'medium', 'charge': 'neutral', 'hydrophobic': False, 'aromatic': False},
    'D': {'name': 'Aspartic acid', 'type': 'acidic', 'size': 'medium', 'charge': 'negative', 'hydrophobic': False, 'aromatic': False},
    'C': {'name': 'Cysteine', 'type': 'polar', 'size': 'small', 'charge': 'neutral', 'hydrophobic': False, 'aromatic': False},
    'Q': {'name': 'Glutamine', 'type': 'polar', 'size': 'medium', 'charge': 'neutral', 'hydrophobic': False, 'aromatic': False},
    'E': {'name': 'Glutamic acid', 'type': 'acidic', 'size': 'medium', 'charge': 'negative', 'hydrophobic': False, 'aromatic': False},
    'G': {'name': 'Glycine', 'type': 'nonpolar', 'size': 'tiny', 'charge': 'neutral', 'hydrophobic': True, 'aromatic': False},
    'H': {'name': 'Histidine', 'type': 'basic', 'size': 'large', 'charge': 'positive', 'hydrophobic': False, 'aromatic': True},
    'I': {'name': 'Isoleucine', 'type': 'nonpolar', 'size': 'medium', 'charge': 'neutral', 'hydrophobic': True, 'aromatic': False},
    'L': {'name': 'Leucine', 'type': 'nonpolar', 'size': 'medium', 'charge': 'neutral', 'hydrophobic': True, 'aromatic': False},
    'K': {'name': 'Lysine', 'type': 'basic', 'size': 'large', 'charge': 'positive', 'hydrophobic': False, 'aromatic': False},
    'M': {'name': 'Methionine', 'type': 'nonpolar', 'size': 'medium', 'charge': 'neutral', 'hydrophobic': True, 'aromatic': False},
    'F': {'name': 'Phenylalanine', 'type': 'nonpolar', 'size': 'large', 'charge': 'neutral', 'hydrophobic': True, 'aromatic': True},
    'P': {'name': 'Proline', 'type': 'nonpolar', 'size': 'small', 'charge': 'neutral', 'hydrophobic': True, 'aromatic': False},
    'S': {'name': 'Serine', 'type': 'polar', 'size': 'small', 'charge': 'neutral', 'hydrophobic': False, 'aromatic': False},
    'T': {'name': 'Threonine', 'type': 'polar', 'size': 'small', 'charge': 'neutral', 'hydrophobic': False, 'aromatic': False},
    'W': {'name': 'Tryptophan', 'type': 'nonpolar', 'size': 'bulky', 'charge': 'neutral', 'hydrophobic': True, 'aromatic': True},
    'Y': {'name': 'Tyrosine', 'type': 'polar', 'size': 'large', 'charge': 'neutral', 'hydrophobic': True, 'aromatic': True},
    'V': {'name': 'Valine', 'type': 'nonpolar', 'size': 'small', 'charge': 'neutral', 'hydrophobic': True, 'aromatic': False},
}

# ============================================================================
# DATA LOADING AND PREPROCESSING
# ============================================================================

def load_pyampa_data(filepath):
    """Load and validate PYAMPA data."""
    global GOOD_AMP_THRESHOLD, GOOD_HEMO_THRESHOLD, GOOD_FITNESS_THRESHOLD
    
    try:
        df = pd.read_csv(filepath)
        
        # Flexible column name matching
        column_mapping = {}
        
        # Required columns with possible variations
        column_variations = {
            'Position': ['Position', 'position', 'Pos'],
            'Original_AA': ['Original_AA', 'original_aa', 'OriginalAA', 'Original'],
            'Mutated_AA': ['Mutated_AA', 'mutated_aa', 'MutatedAA', 'Mutated'],
            'AMP_Prob': ['AMP_Prob', 'AMP_Probab', 'AMP_Probability', 'amp_prob'],
            'Hemo_Prob': ['Hemo_Prob', 'Hemolytic_P', 'Hemolytic_Prob', 'Hemolytic_Probability', 'hemo_prob'],
            'Fitness': ['Fitness', 'Fitness_Score', 'fitness', 'fitness_score']
        }
        
        # Find matching columns
        for standard_name, variations in column_variations.items():
            for var in variations:
                if var in df.columns:
                    column_mapping[var] = standard_name
                    break
        
        # Rename columns to standard names
        df = df.rename(columns=column_mapping)
        
        # Check if all required columns are now present
        required_cols = ['Position', 'Original_AA', 'Mutated_AA', 'AMP_Prob', 'Hemo_Prob', 'Fitness']
        missing = [col for col in required_cols if col not in df.columns]
        
        if missing:
            print(f"❌ Missing columns: {missing}")
            print(f"   Available columns in your CSV: {list(df.columns)}")
            print(f"\n   Please ensure your CSV has these columns:")
            print(f"   - Position (or 'Pos')")
            print(f"   - Original_AA (or 'Original')")
            print(f"   - Mutated_AA (or 'Mutated')")
            print(f"   - AMP_Prob (or 'AMP_Probab' or 'AMP_Probability')")
            print(f"   - Hemo_Prob (or 'Hemolytic_P' or 'Hemolytic_Prob')")
            print(f"   - Fitness (or 'Fitness_Score')")
            return None
        
        print(f"✅ Loaded columns successfully:")
        for old, new in column_mapping.items():
            if old != new:
                print(f"   '{old}' → '{new}'")
        
        # Calculate thresholds
        GOOD_AMP_THRESHOLD = df['AMP_Prob'].median()
        GOOD_HEMO_THRESHOLD = df['Hemo_Prob'].median()
        # GOOD_FITNESS_THRESHOLD already set to 0.0 (biologically meaningful)
        
        print(f"\n📊 Analysis Thresholds:")
        print(f"   AMP Probability    ≥ {GOOD_AMP_THRESHOLD:.4f} (median, 50th percentile)")
        print(f"   Hemolytic Prob     ≤ {GOOD_HEMO_THRESHOLD:.4f} (median, 50th percentile)")
        print(f"   Fitness Score      ≥ {GOOD_FITNESS_THRESHOLD:.4f} (zero baseline)")
        print(f"\n   Note: Median thresholds divide mutations into balanced groups.")
        print(f"         Zero fitness threshold represents wildtype performance.")
        
        # Add property columns
        df['Original_Type'] = df['Original_AA'].map(lambda x: AA_PROPERTIES.get(x, {}).get('type', 'unknown'))
        df['Mutated_Type'] = df['Mutated_AA'].map(lambda x: AA_PROPERTIES.get(x, {}).get('type', 'unknown'))
        df['Original_Size'] = df['Original_AA'].map(lambda x: AA_PROPERTIES.get(x, {}).get('size', 'unknown'))
        df['Mutated_Size'] = df['Mutated_AA'].map(lambda x: AA_PROPERTIES.get(x, {}).get('size', 'unknown'))
        df['Original_Charge'] = df['Original_AA'].map(lambda x: AA_PROPERTIES.get(x, {}).get('charge', 'unknown'))
        df['Mutated_Charge'] = df['Mutated_AA'].map(lambda x: AA_PROPERTIES.get(x, {}).get('charge', 'unknown'))
        df['Original_Hydrophobic'] = df['Original_AA'].map(lambda x: AA_PROPERTIES.get(x, {}).get('hydrophobic', False))
        df['Mutated_Hydrophobic'] = df['Mutated_AA'].map(lambda x: AA_PROPERTIES.get(x, {}).get('hydrophobic', False))
        df['Original_Aromatic'] = df['Original_AA'].map(lambda x: AA_PROPERTIES.get(x, {}).get('aromatic', False))
        df['Mutated_Aromatic'] = df['Mutated_AA'].map(lambda x: AA_PROPERTIES.get(x, {}).get('aromatic', False))
        
        # Calculate property changes
        df['Type_Change'] = df['Original_Type'] + '→' + df['Mutated_Type']
        df['Charge_Change'] = df['Original_Charge'] + '→' + df['Mutated_Charge']
        df['Hydrophobic_Change'] = df.apply(
            lambda r: 'hydrophobic→polar' if r['Original_Hydrophobic'] and not r['Mutated_Hydrophobic']
            else 'polar→hydrophobic' if not r['Original_Hydrophobic'] and r['Mutated_Hydrophobic']
            else 'no_change', axis=1
        )
        df['Aromatic_Change'] = df.apply(
            lambda r: 'aromatic→non-aromatic' if r['Original_Aromatic'] and not r['Mutated_Aromatic']
            else 'non-aromatic→aromatic' if not r['Original_Aromatic'] and r['Mutated_Aromatic']
            else 'no_change', axis=1
        )
        
        # Label mutations as good/bad using calculated thresholds
        df['Good_AMP'] = df['AMP_Prob'] >= GOOD_AMP_THRESHOLD
        df['Good_Hemo'] = df['Hemo_Prob'] <= GOOD_HEMO_THRESHOLD
        df['Good_Fitness'] = df['Fitness'] >= GOOD_FITNESS_THRESHOLD
        df['Overall_Good'] = df['Good_AMP'] & df['Good_Hemo'] & df['Good_Fitness']
        
        return df
        
    except Exception as e:
        print(f"❌ Error loading data: {e}")
        return None


# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

def analyze_amino_acid_preferences(df):
    """Analyze which amino acids are generally beneficial."""
    
    print("\n" + "="*70)
    print("AMINO ACID PREFERENCE ANALYSIS")
    print("="*70)
    
    # Group by mutated amino acid
    aa_stats = df.groupby('Mutated_AA').agg({
        'AMP_Prob': ['mean', 'std', 'count'],
        'Hemo_Prob': ['mean', 'std'],
        'Fitness': ['mean', 'std'],
        'Good_AMP': 'sum',
        'Good_Hemo': 'sum',
        'Overall_Good': 'sum'
    }).round(3)
    
    aa_stats.columns = ['_'.join(col).strip() for col in aa_stats.columns.values]
    aa_stats = aa_stats.sort_values('AMP_Prob_mean', ascending=False)
    
    print("\n📊 Top Amino Acids by Antimicrobial Activity:")
    print(aa_stats[['AMP_Prob_mean', 'Hemo_Prob_mean', 'Fitness_mean', 'AMP_Prob_count']].head(10).to_string())
    
    print("\n📊 Top Amino Acids by Safety (Low Hemolysis):")
    safe_aa = aa_stats.sort_values('Hemo_Prob_mean', ascending=True)
    print(safe_aa[['AMP_Prob_mean', 'Hemo_Prob_mean', 'Fitness_mean', 'AMP_Prob_count']].head(10).to_string())
    
    print("\n📊 Top Amino Acids by Overall Fitness:")
    fit_aa = aa_stats.sort_values('Fitness_mean', ascending=False)
    print(fit_aa[['AMP_Prob_mean', 'Hemo_Prob_mean', 'Fitness_mean', 'AMP_Prob_count']].head(10).to_string())
    
    return aa_stats


def analyze_property_changes(df):
    """Analyze what types of substitutions are beneficial."""
    
    print("\n" + "="*70)
    print("PROPERTY CHANGE ANALYSIS")
    print("="*70)
    
    # Type changes
    print("\n🔬 Amino Acid Type Changes:")
    type_changes = df.groupby('Type_Change').agg({
        'AMP_Prob': 'mean',
        'Hemo_Prob': 'mean',
        'Fitness': 'mean',
        'Overall_Good': 'sum',
        'Position': 'count'
    }).round(3)
    type_changes.columns = ['Avg_AMP', 'Avg_Hemo', 'Avg_Fitness', 'Good_Count', 'Total_Count']
    type_changes = type_changes.sort_values('Avg_Fitness', ascending=False)
    print(type_changes.head(15).to_string())
    
    # Hydrophobicity changes
    print("\n💧 Hydrophobicity Changes:")
    hydro_changes = df[df['Hydrophobic_Change'] != 'no_change'].groupby('Hydrophobic_Change').agg({
        'AMP_Prob': 'mean',
        'Hemo_Prob': 'mean',
        'Fitness': 'mean',
        'Overall_Good': 'sum',
        'Position': 'count'
    }).round(3)
    hydro_changes.columns = ['Avg_AMP', 'Avg_Hemo', 'Avg_Fitness', 'Good_Count', 'Total_Count']
    print(hydro_changes.to_string())
    
    # Aromatic changes
    print("\n🔷 Aromatic Changes:")
    arom_changes = df[df['Aromatic_Change'] != 'no_change'].groupby('Aromatic_Change').agg({
        'AMP_Prob': 'mean',
        'Hemo_Prob': 'mean',
        'Fitness': 'mean',
        'Overall_Good': 'sum',
        'Position': 'count'
    }).round(3)
    arom_changes.columns = ['Avg_AMP', 'Avg_Hemo', 'Avg_Fitness', 'Good_Count', 'Total_Count']
    print(arom_changes.to_string())
    
    return type_changes, hydro_changes, arom_changes


def analyze_position_patterns(df):
    """Analyze position-specific patterns."""
    
    print("\n" + "="*70)
    print("POSITION-SPECIFIC ANALYSIS")
    print("="*70)
    
    positions = sorted(df['Position'].unique())
    
    for pos in positions:
        pos_data = df[df['Position'] == pos]
        
        if len(pos_data) < 3:  # Skip positions with too few mutations
            continue
        
        original_aa = pos_data['Original_AA'].iloc[0]
        
        print(f"\n🔍 Position {pos} (Original: {original_aa} - {AA_PROPERTIES.get(original_aa, {}).get('name', 'Unknown')})")
        print(f"   Type: {AA_PROPERTIES.get(original_aa, {}).get('type', 'unknown')}, "
              f"Charge: {AA_PROPERTIES.get(original_aa, {}).get('charge', 'unknown')}")
        
        # Best mutations at this position
        best_muts = pos_data.nlargest(3, 'Fitness')[['Mutated_AA', 'AMP_Prob', 'Hemo_Prob', 'Fitness']]
        print("\n   Best mutations:")
        for idx, row in best_muts.iterrows():
            mut_aa = row['Mutated_AA']
            print(f"   → {original_aa}→{mut_aa}: AMP={row['AMP_Prob']:.3f}, Hemo={row['Hemo_Prob']:.3f}, Fit={row['Fitness']:+.3f}")
        
        # Identify patterns
        good_muts = pos_data[pos_data['Overall_Good']]
        if len(good_muts) > 0:
            good_types = good_muts['Mutated_Type'].value_counts()
            print(f"\n   ✅ Beneficial mutation types: {', '.join([f'{t}({c})' for t, c in good_types.items()])}")


def extract_principles(df, aa_stats, type_changes, hydro_changes, arom_changes):
    """Extract actionable design principles."""
    
    print("\n" + "="*70)
    print("EXTRACTED DESIGN PRINCIPLES")
    print("="*70)
    
    principles = []
    
    # Principle 1: Best amino acids overall
    best_aa = aa_stats.sort_values('Fitness_mean', ascending=False).head(5)
    principle_1 = {
        "principle": "Preferred Amino Acid Substitutions",
        "rule": f"Generally beneficial amino acids: {', '.join(best_aa.index.tolist())}",
        "rationale": "These amino acids show highest average fitness across all positions",
        "data": best_aa[['AMP_Prob_mean', 'Hemo_Prob_mean', 'Fitness_mean']].to_dict()
    }
    principles.append(principle_1)
    print(f"\n✅ Principle 1: {principle_1['principle']}")
    print(f"   {principle_1['rule']}")
    print(f"   Rationale: {principle_1['rationale']}")
    
    # Principle 2: Hydrophobicity changes
    if len(hydro_changes) > 0:
        best_hydro = hydro_changes.loc[hydro_changes['Avg_Fitness'].idxmax()]
        principle_2 = {
            "principle": "Hydrophobicity Modulation",
            "rule": f"Best change: {hydro_changes['Avg_Fitness'].idxmax()} (Fitness: {best_hydro['Avg_Fitness']:+.3f})",
            "rationale": "Modulating hydrophobicity affects membrane interaction",
            "data": hydro_changes.to_dict()
        }
        principles.append(principle_2)
        print(f"\n✅ Principle 2: {principle_2['principle']}")
        print(f"   {principle_2['rule']}")
        print(f"   Rationale: {principle_2['rationale']}")
    
    # Principle 3: Aromatic changes
    if len(arom_changes) > 0:
        best_arom = arom_changes.loc[arom_changes['Avg_Fitness'].idxmax()]
        principle_3 = {
            "principle": "Aromatic Residue Effects",
            "rule": f"Best change: {arom_changes['Avg_Fitness'].idxmax()} (Fitness: {best_arom['Avg_Fitness']:+.3f})",
            "rationale": "Aromatic residues affect membrane insertion and stability",
            "data": arom_changes.to_dict()
        }
        principles.append(principle_3)
        print(f"\n✅ Principle 3: {principle_3['principle']}")
        print(f"   {principle_3['rule']}")
        print(f"   Rationale: {principle_3['rationale']}")
    
    # Principle 4: Type substitutions
    best_type_changes = type_changes.nlargest(5, 'Avg_Fitness')
    principle_4 = {
        "principle": "Beneficial Amino Acid Type Substitutions",
        "rule": f"Top changes: {', '.join(best_type_changes.index.tolist())}",
        "rationale": "These substitution patterns consistently improve peptide properties",
        "data": best_type_changes[['Avg_AMP', 'Avg_Hemo', 'Avg_Fitness']].to_dict()
    }
    principles.append(principle_4)
    print(f"\n✅ Principle 4: {principle_4['principle']}")
    print(f"   {principle_4['rule']}")
    
    # Principle 5: Safety-focused
    safe_aa = aa_stats.nsmallest(5, 'Hemo_Prob_mean')
    principle_5 = {
        "principle": "Safety-Focused Substitutions",
        "rule": f"Lowest hemolysis: {', '.join(safe_aa.index.tolist())}",
        "rationale": "These amino acids minimize hemolytic activity",
        "data": safe_aa[['AMP_Prob_mean', 'Hemo_Prob_mean', 'Fitness_mean']].to_dict()
    }
    principles.append(principle_5)
    print(f"\n✅ Principle 5: {principle_5['principle']}")
    print(f"   {principle_5['rule']}")
    print(f"   Rationale: {principle_5['rationale']}")
    
    # Principle 6: Activity-focused
    active_aa = aa_stats.nlargest(5, 'AMP_Prob_mean')
    principle_6 = {
        "principle": "Activity-Focused Substitutions",
        "rule": f"Highest antimicrobial activity: {', '.join(active_aa.index.tolist())}",
        "rationale": "These amino acids maximize antimicrobial potency",
        "data": active_aa[['AMP_Prob_mean', 'Hemo_Prob_mean', 'Fitness_mean']].to_dict()
    }
    principles.append(principle_6)
    print(f"\n✅ Principle 6: {principle_6['principle']}")
    print(f"   {principle_6['rule']}")
    print(f"   Rationale: {principle_6['rationale']}")
    
    return principles


def generate_application_guide(df, principles):
    """Generate practical application guide."""
    
    print("\n" + "="*70)
    print("PRACTICAL APPLICATION GUIDE")
    print("="*70)
    
    print("\n📖 How to Apply These Principles to New Sequences:")
    print("\n1. IDENTIFY SEQUENCE CONTEXT")
    print("   - Determine amino acid properties at target position")
    print("   - Consider surrounding residues")
    
    print("\n2. SELECT MUTATION STRATEGY")
    print("   - For HIGH ACTIVITY: Use amino acids from Principle 6")
    print("   - For LOW HEMOLYSIS: Use amino acids from Principle 5")
    print("   - For BALANCE: Use amino acids from Principle 1")
    
    print("\n3. APPLY PROPERTY-BASED RULES")
    
    # Hydrophobicity rule
    hydro_change = df.groupby('Hydrophobic_Change').agg({'Fitness': 'mean'}).round(3)
    if 'hydrophobic→polar' in hydro_change.index and 'polar→hydrophobic' in hydro_change.index:
        if hydro_change.loc['hydrophobic→polar', 'Fitness'] > hydro_change.loc['polar→hydrophobic', 'Fitness']:
            print("   - HYDROPHOBICITY: Replacing hydrophobic with polar residues is generally beneficial")
        else:
            print("   - HYDROPHOBICITY: Replacing polar with hydrophobic residues is generally beneficial")
    
    # Aromatic rule
    arom_change = df.groupby('Aromatic_Change').agg({'Fitness': 'mean'}).round(3)
    if 'aromatic→non-aromatic' in arom_change.index and 'non-aromatic→aromatic' in arom_change.index:
        if arom_change.loc['aromatic→non-aromatic', 'Fitness'] > arom_change.loc['non-aromatic→aromatic', 'Fitness']:
            print("   - AROMATIC: Replacing aromatic with non-aromatic residues is generally beneficial")
        else:
            print("   - AROMATIC: Replacing non-aromatic with aromatic residues is generally beneficial")
    
    print("\n4. EXAMPLE APPLICATION")
    print("   Given sequence: KRRFHWWWMFLRR")
    print("   Target: Reduce hemolysis at position 7 (W)")
    print("")
    print("   Step 1: Identify W properties → aromatic, hydrophobic, bulky")
    print("   Step 2: Apply Principle 5 (safety) → Try G, A, or S")
    print("   Step 3: Apply hydrophobicity rule → hydrophobic→polar beneficial")
    print("   Step 4: Select G (small, polar, non-aromatic)")
    print("   Result: KRRFHWGWMFLRR (predicted lower hemolysis)")


def save_results(principles, aa_stats, output_file):
    """Save principles to JSON file."""
    
    # Convert DataFrames to dictionaries for JSON serialization
    output_data = {
        "principles": principles,
        "amino_acid_rankings": {
            "by_fitness": aa_stats.sort_values('Fitness_mean', ascending=False).head(10).to_dict(),
            "by_amp": aa_stats.sort_values('AMP_Prob_mean', ascending=False).head(10).to_dict(),
            "by_safety": aa_stats.sort_values('Hemo_Prob_mean', ascending=True).head(10).to_dict(),
        },
        "thresholds": {
            "good_amp": GOOD_AMP_THRESHOLD,
            "good_hemo": GOOD_HEMO_THRESHOLD,
            "good_fitness": GOOD_FITNESS_THRESHOLD,
            "threshold_method": {
                "amp": "median",
                "hemo": "median", 
                "fitness": "zero_baseline"
            }
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n✅ Principles saved to: {output_file}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main analysis workflow."""
    
    print("\n" + "="*70)
    print("PYAMPA PATTERN ANALYZER")
    print("="*70)
    print("Extracting Structure-Activity Relationship (SAR) Principles")
    print("="*70)
    
    # Load data
    print("\nLoading PYAMPA data...")
    df = load_pyampa_data(PYAMPA_FILE)
    
    if df is None:
        print("\n❌ Failed to load data. Please check your CSV file.")
        return
    
    print(f"✅ Loaded {len(df)} mutations across {df['Position'].nunique()} positions")
    print(f"✅ Analyzing {df['Original_AA'].nunique()} → {df['Mutated_AA'].nunique()} amino acid changes")
    
    # Run analyses
    aa_stats = analyze_amino_acid_preferences(df)
    type_changes, hydro_changes, arom_changes = analyze_property_changes(df)
    analyze_position_patterns(df)
    principles = extract_principles(df, aa_stats, type_changes, hydro_changes, arom_changes)
    generate_application_guide(df, principles)
    
    # Save results
    save_results(principles, aa_stats, OUTPUT_PRINCIPLES)
    
    # Summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\n📊 Summary Statistics:")
    print(f"   Total mutations analyzed: {len(df)}")
    print(f"   Beneficial mutations (all criteria): {df['Overall_Good'].sum()}")
    print(f"   High AMP mutations: {df['Good_AMP'].sum()}")
    print(f"   Low hemolysis mutations: {df['Good_Hemo'].sum()}")
    print(f"   Positive fitness mutations: {df['Good_Fitness'].sum()}")
    print(f"\n📋 Extracted {len(principles)} design principles")
    print(f"✅ Results saved to: {OUTPUT_PRINCIPLES}")
    print("\n" + "="*70 + "\n")


if __name__ == "__main__":
    main()
