"""
Smart Peptide Variant Generator
================================
Generates optimized peptide variants based on PYAMPA-derived design principles.

Features:
- Specify which positions to modify (e.g., hemolytic hotspots from HemoPI2)
- Apply data-driven mutation strategies
- Generate multiple candidates with different mutation combinations
- Predict HC50 and MIC changes
- Rank by safety improvement

"""

import pandas as pd
import numpy as np
import json
from itertools import combinations
import sys

# ============================================================================
# CONFIGURATION
# ============================================================================

PRINCIPLES_FILE = "principles/design_principles_YOURSEQ.json"
OUTPUT_FILE = "YOURSEQ_variants.csv"

# Prediction parameters
#BASE_HC50 = number(3 decimal)  # Default baseline HC50 (can be changed)
#BASE_MIC = number(3 decimal)   # Default baseline MIC (can be changed)

# Generation parameters
MAX_MUTATIONS = 2        # Maximum simultaneous mutations
TOP_N_CANDIDATES = 10    # Number of candidates to generate

# Optimization strategy
#STRATEGY = "safety"      # Options: "safety", "balanced", "activity"

# ============================================================================
# AMINO ACID PROPERTY DATABASE
# ============================================================================

AA_PROPERTIES = {
    'A': {'name': 'Alanine', 'type': 'nonpolar', 'hydrophobic': True},
    'R': {'name': 'Arginine', 'type': 'basic', 'hydrophobic': False},
    'N': {'name': 'Asparagine', 'type': 'polar', 'hydrophobic': False},
    'D': {'name': 'Aspartic acid', 'type': 'acidic', 'hydrophobic': False},
    'C': {'name': 'Cysteine', 'type': 'polar', 'hydrophobic': False},
    'Q': {'name': 'Glutamine', 'type': 'polar', 'hydrophobic': False},
    'E': {'name': 'Glutamic acid', 'type': 'acidic', 'hydrophobic': False},
    'G': {'name': 'Glycine', 'type': 'nonpolar', 'hydrophobic': False},
    'H': {'name': 'Histidine', 'type': 'basic', 'hydrophobic': False},
    'I': {'name': 'Isoleucine', 'type': 'nonpolar', 'hydrophobic': True},
    'L': {'name': 'Leucine', 'type': 'nonpolar', 'hydrophobic': True},
    'K': {'name': 'Lysine', 'type': 'basic', 'hydrophobic': False},
    'M': {'name': 'Methionine', 'type': 'nonpolar', 'hydrophobic': True},
    'F': {'name': 'Phenylalanine', 'type': 'nonpolar', 'hydrophobic': True},
    'P': {'name': 'Proline', 'type': 'nonpolar', 'hydrophobic': True},
    'S': {'name': 'Serine', 'type': 'polar', 'hydrophobic': False},
    'T': {'name': 'Threonine', 'type': 'polar', 'hydrophobic': False},
    'W': {'name': 'Tryptophan', 'type': 'nonpolar', 'hydrophobic': True},
    'Y': {'name': 'Tyrosine', 'type': 'polar', 'hydrophobic': True},
    'V': {'name': 'Valine', 'type': 'nonpolar', 'hydrophobic': True},
}

# ============================================================================
# LOAD DESIGN PRINCIPLES
# ============================================================================

def load_principles(filepath):
    """Load design principles from JSON file."""
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
        print(f"✅ Loaded design principles from {filepath}")
        return data
    except FileNotFoundError:
        print(f"❌ Error: {filepath} not found!")
        print(f"   Please run pyampa_pattern_analyzer.py first to generate principles.")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error loading principles: {e}")
        sys.exit(1)


def get_safety_ranked_amino_acids(principles):
    """Get amino acids ranked by safety (low hemolysis)."""
    safety_data = principles['amino_acid_rankings']['by_safety']
    aa_list = list(safety_data['Hemo_Prob_mean'].keys())
    hemo_values = safety_data['Hemo_Prob_mean']
    
    # Sort by hemolysis (lower is better)
    sorted_aa = sorted(aa_list, key=lambda x: hemo_values[x])
    return sorted_aa, safety_data


def get_activity_ranked_amino_acids(principles):
    """Get amino acids ranked by activity (high AMP)."""
    activity_data = principles['amino_acid_rankings']['by_amp']
    aa_list = list(activity_data['AMP_Prob_mean'].keys())
    amp_values = activity_data['AMP_Prob_mean']
    
    # Sort by AMP (higher is better)
    sorted_aa = sorted(aa_list, key=lambda x: amp_values[x], reverse=True)
    return sorted_aa, activity_data


def get_fitness_ranked_amino_acids(principles):
    """Get amino acids ranked by fitness."""
    fitness_data = principles['amino_acid_rankings']['by_fitness']
    aa_list = list(fitness_data['Fitness_mean'].keys())
    fitness_values = fitness_data['Fitness_mean']
    
    # Sort by fitness (higher is better)
    sorted_aa = sorted(aa_list, key=lambda x: fitness_values[x], reverse=True)
    return sorted_aa, fitness_data


# ============================================================================
# MUTATION GENERATION
# ============================================================================

def apply_mutation(sequence, position, amino_acid):
    """Apply a single mutation to a sequence."""
    seq_list = list(sequence)
    seq_list[position - 1] = amino_acid
    return ''.join(seq_list)


def get_mutation_candidates(sequence, positions, principles, strategy="safety", top_n=5):
    """
    Get top candidate amino acids for each position based on strategy.
    
    Args:
        sequence: Original peptide sequence
        positions: List of positions to mutate (1-indexed)
        principles: Loaded design principles
        strategy: "safety", "balanced", or "activity"
        top_n: Number of candidate amino acids per position
        
    Returns:
        Dictionary: {position: [list of candidate amino acids]}
    """
    
    if strategy == "safety":
        ranked_aa, aa_data = get_safety_ranked_amino_acids(principles)
        print(f"\n📋 Using SAFETY strategy (minimize hemolysis)")
        print(f"   Prioritized amino acids: {', '.join(ranked_aa[:5])}")
        
    elif strategy == "activity":
        ranked_aa, aa_data = get_activity_ranked_amino_acids(principles)
        print(f"\n📋 Using ACTIVITY strategy (maximize AMP)")
        print(f"   Prioritized amino acids: {', '.join(ranked_aa[:5])}")
        
    else:  # balanced
        ranked_aa, aa_data = get_fitness_ranked_amino_acids(principles)
        print(f"\n📋 Using BALANCED strategy (optimize overall fitness)")
        print(f"   Prioritized amino acids: {', '.join(ranked_aa[:5])}")
    
    candidates = {}
    
    for pos in positions:
        original_aa = sequence[pos - 1]
        
        # Get top amino acids that are different from original
        position_candidates = [aa for aa in ranked_aa if aa != original_aa][:top_n]
        
        candidates[pos] = position_candidates
        
        print(f"\n   Position {pos} ({original_aa}):")
        print(f"   Candidates: {', '.join(position_candidates)}")
    
    return candidates


def generate_variant_combinations(sequence, position_candidates, max_mutations=4):
    """
    Generate all possible mutation combinations.
    
    Args:
        sequence: Original sequence
        position_candidates: Dict of {position: [candidate AAs]}
        max_mutations: Maximum number of simultaneous mutations
        
    Returns:
        List of variant dictionaries
    """
    
    variants = []
    positions = list(position_candidates.keys())
    
    # Generate single mutations
    for pos in positions:
        for aa in position_candidates[pos]:
            mutation_dict = {pos: aa}
            new_seq = sequence
            for p, a in mutation_dict.items():
                new_seq = apply_mutation(new_seq, p, a)
            
            variants.append({
                'sequence': new_seq,
                'mutations': mutation_dict,
                'mutation_string': f"{sequence[pos-1]}{pos}{aa}",
                'num_mutations': 1
            })
    
    # Generate double mutations
    if max_mutations >= 2 and len(positions) >= 2:
        for pos1, pos2 in combinations(positions, 2):
            # Take top 2 candidates per position for combinations
            for aa1 in position_candidates[pos1][:2]:
                for aa2 in position_candidates[pos2][:2]:
                    mutation_dict = {pos1: aa1, pos2: aa2}
                    new_seq = sequence
                    for p, a in mutation_dict.items():
                        new_seq = apply_mutation(new_seq, p, a)
                    
                    mut_str = f"{sequence[pos1-1]}{pos1}{aa1}+{sequence[pos2-1]}{pos2}{aa2}"
                    
                    variants.append({
                        'sequence': new_seq,
                        'mutations': mutation_dict,
                        'mutation_string': mut_str,
                        'num_mutations': 2
                    })
    
    # Generate triple mutations
    if max_mutations >= 3 and len(positions) >= 3:
        for pos1, pos2, pos3 in combinations(positions, 3):
            # Take top 1 candidate per position for triple combinations
            aa1 = position_candidates[pos1][0]
            aa2 = position_candidates[pos2][0]
            aa3 = position_candidates[pos3][0]
            
            mutation_dict = {pos1: aa1, pos2: aa2, pos3: aa3}
            new_seq = sequence
            for p, a in mutation_dict.items():
                new_seq = apply_mutation(new_seq, p, a)
            
            mut_str = f"{sequence[pos1-1]}{pos1}{aa1}+{sequence[pos2-1]}{pos2}{aa2}+{sequence[pos3-1]}{pos3}{aa3}"
            
            variants.append({
                'sequence': new_seq,
                'mutations': mutation_dict,
                'mutation_string': mut_str,
                'num_mutations': 3
            })

    if max_mutations >= 4 and len(positions) >= 4:
        for combo in combinations(positions, 4):
            pos1, pos2, pos3, pos4 = combo
            for aa1 in position_candidates[pos1]:
                for aa2 in position_candidates[pos2]:
                    for aa3 in position_candidates[pos3]:
                        for aa4 in position_candidates[pos4]:
                            mutation_dict = {pos1: aa1, pos2: aa2, pos3: aa3, pos4: aa4}
                            new_seq = sequence
                            new_seq = apply_mutation(new_seq, pos1, aa1)
                            new_seq = apply_mutation(new_seq, pos2, aa2)
                            new_seq = apply_mutation(new_seq, pos3, aa3)
                            new_seq = apply_mutation(new_seq, pos4, aa4)
                            mut_str = (
                                f"{sequence[pos1-1]}{pos1}{aa1}+"
                                f"{sequence[pos2-1]}{pos2}{aa2}+"
                                f"{sequence[pos3-1]}{pos3}{aa3}+"
                                f"{sequence[pos4-1]}{pos4}{aa4}"
                            )
                            variants.append({
                                'sequence': new_seq,
                                'mutations': mutation_dict,
                                'mutation_string': mut_str,
                                'num_mutations': 4
                            })
    return variants

def generate_single_mutations_only(sequence, position_candidates):
    """Generate ONLY single-mutation variants (to unmask them)."""
    variants = []
    positions = list(position_candidates.keys())
    
    for pos in positions:
        for aa in position_candidates[pos]:
            mutation_dict = {pos: aa}
            new_seq = apply_mutation(sequence, pos, aa)
            
            variants.append({
                'sequence': new_seq,
                'mutations': mutation_dict,
                'mutation_string': f"{sequence[pos-1]}{pos}{aa}",
                'num_mutations': 1
            })
    return variants

# ============================================================================
# PROPERTY PREDICTION
# ============================================================================

def predict_properties(original_seq, mutations, principles, base_hc50, base_mic):
    """
    Predict HC50 and MIC for mutated sequence.
    
    Uses amino acid property data from principles to estimate changes.
    """
    
    safety_data = principles['amino_acid_rankings']['by_safety']
    activity_data = principles['amino_acid_rankings']['by_amp']
    
    hemo_modifier = 1.0
    amp_modifier = 1.0
    
    for pos, new_aa in mutations.items():
        original_aa = original_seq[pos - 1]
        
        # Get hemolysis probabilities
        original_hemo = safety_data['Hemo_Prob_mean'].get(original_aa, 0.90)
        new_hemo = safety_data['Hemo_Prob_mean'].get(new_aa, 0.90)
        
        # Get AMP probabilities
        original_amp = activity_data['AMP_Prob_mean'].get(original_aa, 0.75)
        new_amp = activity_data['AMP_Prob_mean'].get(new_aa, 0.75)
        
        # Calculate modifiers
        hemo_modifier *= (new_hemo / original_hemo)
        amp_modifier *= (new_amp / original_amp)
    
    # Predict new values
    # HC50: Lower Hemo_Prob = Higher HC50 (safer), so we DIVIDE
    predicted_hc50 = base_hc50 / hemo_modifier if hemo_modifier > 0 else base_hc50
    
    # MIC: Higher AMP_Prob = Lower MIC (better), so we DIVIDE
    predicted_mic = base_mic / amp_modifier if amp_modifier > 0 else base_mic
    
    # Calculate improvements
    hc50_improvement = ((predicted_hc50 - base_hc50) / base_hc50) * 100
    mic_improvement = ((base_mic - predicted_mic) / base_mic) * 100
    
    return {
        'predicted_hc50': predicted_hc50,
        'predicted_mic': predicted_mic,
        'hc50_improvement': hc50_improvement,
        'mic_improvement': mic_improvement,
        'hemo_modifier': hemo_modifier,
        'amp_modifier': amp_modifier
    }


# ============================================================================
# RANKING AND SCORING
# ============================================================================

def score_variants(variants, strategy="safety"):
    """
    Score variants based on optimization strategy.
    
    Args:
        variants: List of variant dictionaries with predictions
        strategy: "safety", "balanced", or "activity"
        
    Returns:
        Sorted list of variants with scores
    """
    
    for variant in variants:
        if strategy == "safety":
            # Prioritize HC50 improvement heavily
            score = (variant['hc50_improvement'] * 3.0 + 
                    variant['mic_improvement'] * 0.5)
            
        elif strategy == "activity":
            # Prioritize MIC improvement heavily
            score = (variant['mic_improvement'] * 3.0 + 
                    variant['hc50_improvement'] * 0.5)
            
        else:  # balanced
            # Equal weight to both
            score = (variant['hc50_improvement'] * 1.5 + 
                    variant['mic_improvement'] * 1.5)
        
        # Penalty for too many mutations (prefer simpler)
        score -= (variant['num_mutations'] - 1) * 5
        
        variant['score'] = score
    
    # Sort by score (highest first)
    variants_sorted = sorted(variants, key=lambda x: x['score'], reverse=True)
    
    return variants_sorted


# ============================================================================
# OUTPUT AND DISPLAY
# ============================================================================

def display_results(variants, original_seq, base_hc50, base_mic):
    """Display top variant recommendations."""
    
    print("\n" + "="*80)
    print(f"TOP {len(variants)} RECOMMENDED VARIANTS")
    print("="*80)
    
    print(f"\nOriginal Sequence: {original_seq}")
    print(f"Baseline HC50: {base_hc50:.1f} μM")
    print(f"Baseline MIC:  {base_mic:.3f} μg/mL")
    print()
    
    for idx, variant in enumerate(variants, 1):
        print(f"\n{'─'*80}")
        print(f"Rank #{idx}")
        print(f"{'─'*80}")
        print(f"Sequence:     {variant['sequence']}")
        print(f"Mutations:    {variant['mutation_string']}")
        print(f"")
        print(f"Predicted HC50:  {variant['predicted_hc50']:.1f} μM ({variant['hc50_improvement']:+.1f}%)")
        print(f"Predicted MIC:   {variant['predicted_mic']:.3f} μg/mL ({variant['mic_improvement']:+.1f}%)")
        print(f"Score:           {variant['score']:.2f}")
        
        # Status indicators
        status = []
        if variant['hc50_improvement'] > 0:
            status.append("✅ Improved safety")
        if variant['mic_improvement'] > 0:
            status.append("✅ Improved activity")
        if variant['predicted_hc50'] > 50:
            status.append("⭐ Non-hemolytic predicted")
        
        if status:
            print(f"Status:          {' | '.join(status)}")


def export_results(variants, original_seq, output_file):
    """Export results to CSV."""
    
    export_data = []
    for idx, variant in enumerate(variants, 1):
        export_data.append({
            'rank': idx,
            'sequence': variant['sequence'],
            'original': original_seq,
            'mutations': variant['mutation_string'],
            'num_mutations': variant['num_mutations'],
            'predicted_hc50': round(variant['predicted_hc50'], 2),
            'predicted_mic': round(variant['predicted_mic'], 3),
            'hc50_improvement_%': round(variant['hc50_improvement'], 1),
            'mic_improvement_%': round(variant['mic_improvement'], 1),
            'score': round(variant['score'], 2),
            'hemo_modifier': round(variant['hemo_modifier'], 3),
            'amp_modifier': round(variant['amp_modifier'], 3)
        })
    
    df = pd.DataFrame(export_data)
    df.to_csv(output_file, index=False)
    print(f"\n✅ Results exported to: {output_file}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function."""
    
    print("\n" + "="*80)
    print("SMART PEPTIDE VARIANT GENERATOR")
    print("="*80)
    print("Generate optimized variants based on PYAMPA design principles")
    print("="*80)
    
    # Load principles
    principles = load_principles(PRINCIPLES_FILE)
    
    # Get user input
    print("\n" + "─"*80)
    print("INPUT PARAMETERS")
    print("─"*80)
    
    original_sequence = input("\nEnter your peptide sequence: ").strip().upper()
    if not original_sequence:
        print("❌ No sequence provided!")
        sys.exit(1)
    
    print(f"\n✅ Input sequence: {original_sequence} (Length: {len(original_sequence)})")
    
    # Get baseline values
    print("\n" + "─"*40)
    print("BASELINE VALUES (required every run)")

    while True:
        try:
            base_hc50 = float(input("Baseline HC50 (μM) [higher = safer, e.g. 37.146]: ").strip())
            break
        except ValueError:
            print("   → Please enter a number (e.g. 42.5)") 
    while True:
        try:
            base_mic = float(input("Baseline MIC (μg/mL) [lower = better, e.g. 0.714]: ").strip())
            break
        except ValueError:
            print("   → Please enter a number (e.g. 1.25)")
    print(f"→ Using: HC50 = {base_hc50:.3f} μM  |  MIC = {base_mic:.3f} μg/mL")
    print("─"*40 + "\n")
    

    # Get positions to mutate
    print("\n" + "─"*40)
    print("Which positions would you like to mutate?")
    print("(e.g., '4,5,7' or '1-3,10' for ranges)")
    positions_input = input("Positions: ").strip()
    
    # Parse positions
    positions = []
    for part in positions_input.split(','):
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            positions.extend(range(start, end + 1))
        else:
            positions.append(int(part))
    
    positions = sorted(list(set(positions)))  # Remove duplicates and sort
    
    # Validate positions
    if not all(1 <= p <= len(original_sequence) for p in positions):
        print("❌ Invalid positions! Must be between 1 and sequence length.")
        sys.exit(1)
    
    print(f"✅ Target positions: {positions}")
    print(f"   Original amino acids: {', '.join([f'{p}:{original_sequence[p-1]}' for p in positions])}")
    
    # Get strategy
    print("\n" + "─"*40)
    print("Select optimization strategy:")
    print("  1. Safety (minimize hemolysis) - RECOMMENDED for your case")
    print("  2. Balanced (optimize both)")
    print("  3. Activity (maximize antimicrobial)")
    strategy_choice = input("Choice [1/2/3]: ").strip()
    
    strategy_map = {'1': 'safety', '2': 'balanced', '3': 'activity'}
    strategy = strategy_map.get(strategy_choice, 'safety')
    
    print(f"✅ Selected strategy: {strategy.upper()}")
    

        # Generate candidates
    print("\n" + "="*80)
    print("GENERATING VARIANTS")
    print("="*80)
    
    position_candidates = get_mutation_candidates(
        original_sequence, 
        positions, 
        principles, 
        strategy=strategy,
        top_n=5
    )
    
    # 1. Normal max=2 run (exactly like your original script)
    variants_max2 = generate_variant_combinations(
        original_sequence,
        position_candidates,
        max_mutations=MAX_MUTATIONS
    )
    
    # 2. Single-mutation only (unmasked)
    variants_single = generate_single_mutations_only(
        original_sequence,
        position_candidates
    )
    
    print(f"✅ Generated {len(variants_max2)} max={MAX_MUTATIONS} variants + {len(variants_single)} single-mutation variants")
    
    # Predict properties for both sets
    print(f"\n⏳ Predicting properties...")
    for variant in variants_max2 + variants_single:
        predictions = predict_properties(
            original_sequence,
            variant['mutations'],
            principles,
            base_hc50,
            base_mic
        )
        variant.update(predictions)
    
    # Score and rank each set SEPARATELY (no clash)
    top_max2 = score_variants(variants_max2, strategy=strategy)[:TOP_N_CANDIDATES]
    top_single = score_variants(variants_single, strategy=strategy)[:TOP_N_CANDIDATES]
    
    # Display both clearly
    print("\n" + "="*60)
    print(f"MAX={MAX_MUTATIONS} VARIANTS")
    print("="*60)
    display_results(top_max2, original_sequence, base_hc50, base_mic)
    
    print("\n" + "="*60)
    print("SINGLE-MUTATION VARIANTS (Unmasked)")
    print("="*60)
    display_results(top_single, original_sequence, base_hc50, base_mic)
    
    # Export both sets into ONE CSV file
    output_file = f"{original_sequence}_dual_variants.csv"
    
    # Combine for export (with clear separation)
    combined = []
    for v in top_max2:
        v_copy = v.copy()
        v_copy['group'] = f"max{MAX_MUTATIONS}"
        combined.append(v_copy)
    
    for v in top_single:
        v_copy = v.copy()
        v_copy['group'] = "single"
        combined.append(v_copy)
    
    export_results(combined, original_sequence, output_file)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\n✅ Combined results saved to: {output_file}")
    print(f"   → {len(top_max2)} from max={MAX_MUTATIONS} + {len(top_single)} from single-mutation")
    print("\nNext steps:")
    print("  1. Review both groups in the CSV")
    print("  2. Validate with HemoPI2, ToxinPred, and other tools")
    print("  3. Select 2-3 best candidates for experimental testing")
    print("\n" + "="*80 + "\n")


if __name__ == "__main__":
    main()
