import numpy as np
import hicstraw
import matplotlib.pyplot as plt
import seaborn as sns

# Set resolution
resolution = 500  # Change this to 500 or 1000 as needed

# Define all regions of interest
all_regions = [
    ("E1", 21965574, 21965827),
    ("P16", 21974644, 21975381),
    ("P14_ANRIL", 21994286, 21994909),
    ("P15", 22008685, 22010033),
    ("E2", 22054180, 22054370),
    ("E3", 22097733, 22098066),
    ("E4", 22111884, 22112211),
    ("E5", 22114639, 22115733),
    ("E6", 22117396, 22118792),
    ("E7", 22148199, 22148548),
    ("E8", 22170179, 22170509),
    ("E9", 22211647, 22212006),
    ("MTAP", 21802636, 21937651),
]

# Chromosome
chr_name = "9"

# Define the full region to query (MTAP to furthest downstream)
query_start = (21802636 // resolution) * resolution
query_end = ((max([end for _, _, end in all_regions]) + 10000) // resolution) * resolution

# Function to perform virtual 4C for a given viewpoint
def virtual_4c_for_viewpoint(viewpoint_name, viewpoint_start, viewpoint_end, chr_name, resolution, query_start, query_end):
    """
    Perform virtual 4C analysis for a specific viewpoint
    """
    print(f"\n{'='*80}")
    print(f"Processing viewpoint: {viewpoint_name}")
    print(f"{'='*80}")
    
    # Calculate viewpoint center
    viewpoint_center = (viewpoint_start + viewpoint_end) // 2
    
    # Round coordinates to resolution boundaries
    viewpoint_start_rounded = (viewpoint_start // resolution) * resolution
    viewpoint_end_rounded = ((viewpoint_end // resolution) + 1) * resolution
    
    print(f"Viewpoint: chr{chr_name}:{viewpoint_start}-{viewpoint_end}")
    print(f"Rounded viewpoint: chr{chr_name}:{viewpoint_start_rounded}-{viewpoint_end_rounded}")
    print(f"Full query region: chr{chr_name}:{query_start}-{query_end}")
    
    # Query Hi-C data in two parts (upstream and downstream)
    positions = []
    counts = []
    exclude_distance = 0  # Exclude contacts very close to viewpoint
    
    # Part 1: Upstream (query_start to viewpoint)
    print("Querying upstream region...")
    try:
        result_upstream = hicstraw.straw(
            "observed", 
            'NONE', 
            'data/merged.hic', 
            f'chr{chr_name}:{query_start}:{viewpoint_start_rounded}',
            f'chr{chr_name}:{viewpoint_start_rounded}:{viewpoint_end_rounded}',
            'BP', 
            resolution
        )
        
        # Process upstream results
        for contact in result_upstream:
            if abs(contact.binX - viewpoint_center) > exclude_distance:
                positions.append(contact.binX)
                counts.append(contact.counts)
    except Exception as e:
        print(f"Warning: Error querying upstream region: {e}")
    
    # Part 2: Downstream (viewpoint to query_end)
    print("Querying downstream region...")
    try:
        result_downstream = hicstraw.straw(
            "observed", 
            'NONE', 
            'data/merged.hic', 
            f'chr{chr_name}:{viewpoint_start_rounded}:{viewpoint_end_rounded}',
            f'chr{chr_name}:{viewpoint_start_rounded}:{query_end}',
            'BP', 
            resolution
        )
        
        # Process downstream results
        for contact in result_downstream:
            if abs(contact.binY - viewpoint_center) > exclude_distance:
                positions.append(contact.binY)
                counts.append(contact.counts)
    except Exception as e:
        print(f"Warning: Error querying downstream region: {e}")
    
    # Convert to numpy arrays
    positions = np.array(positions)
    counts = np.array(counts)
    
    print(f"Total contacts retrieved: {len(counts)}")
    if len(counts) > 0:
        print(f"Mean contact frequency: {np.mean(counts):.2f}")
        print(f"Max contact frequency: {np.max(counts):.2f}")
        print(f"Position with max contacts: {positions[np.argmax(counts)]}")
    
    return positions, counts, viewpoint_center, viewpoint_start, viewpoint_end


# Function to create plot
def plot_virtual_4c(positions, counts, viewpoint_name, viewpoint_center, viewpoint_start, viewpoint_end, 
                    chr_name, resolution, all_regions):
    """
    Create virtual 4C plot
    """
    fig, ax = plt.subplots(figsize=(16, 5))
    
    if len(positions) > 0:
        ax.bar(positions, counts, width=resolution, color='steelblue', alpha=0.7, edgecolor='none')
        ax.set_xlabel(f'Position on chr{chr_name} (bp)', fontsize=12)
        ax.set_ylabel('Contact Frequency', fontsize=12)
        ax.set_title(f'Virtual 4C from {viewpoint_name} (chr{chr_name}:{viewpoint_start}-{viewpoint_end})', 
                     fontsize=14, fontweight='bold')
        
        # Add shaded boxes and labels for all regions
        y_max = ax.get_ylim()[1]
        for label, start, end in all_regions:
            # Use red for the current viewpoint, gray for others
            if label == viewpoint_name:
                color = 'red'
                alpha = 0.3
                text_color = 'red'
                text_weight = 'bold'
            else:
                color = 'gray'
                alpha = 0.2
                text_color = 'darkgray'
                text_weight = 'normal'
            
            ax.axvspan(start, end, color=color, alpha=alpha, zorder=0)
            center = (start + end) // 2
            ax.text(center, y_max * 0.9, label, 
                    rotation=90, ha='center', va='top', 
                    fontsize=8, color=text_color, fontweight=text_weight)
        
        # Add viewpoint marker
        ax.axvline(viewpoint_center, color='red', linestyle='-', linewidth=1.5, alpha=0.6, zorder=0)
        
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.ticklabel_format(style='plain', axis='x')
        plt.xticks(rotation=45, ha='right')
    else:
        ax.text(0.5, 0.5, 'No data retrieved', 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
    
    plt.tight_layout()
    
    # Save plot
    plot_filename = f"virtual_4C_{viewpoint_name}_chr{chr_name}_res{resolution}.png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved: {plot_filename}")
    plt.close()


# Function to write WIG file
def write_wig_file(positions, counts, viewpoint_name, viewpoint_start, viewpoint_end, 
                   chr_name, resolution):
    """
    Write virtual 4C data to WIG file
    """
    wig_filename = f"virtual_4C_{viewpoint_name}_chr{chr_name}_res{resolution}.wig"
    print(f"Writing WIG file: {wig_filename}")
    
    if len(positions) > 0:
        # Sort positions and counts together
        sorted_indices = np.argsort(positions)
        sorted_positions = positions[sorted_indices]
        sorted_counts = counts[sorted_indices]
        
        with open(wig_filename, 'w') as f:
            # Write WIG header with viewpoint-specific information
            f.write(f'track type=wiggle_0 name="Virtual_4C_{viewpoint_name}" ')
            f.write(f'description="Virtual 4C from {viewpoint_name} (chr{chr_name}:{viewpoint_start}-{viewpoint_end}) at {resolution}bp resolution"\n')
            
            # Use variableStep format
            f.write(f"variableStep chrom=chr{chr_name} span={resolution}\n")
            
            # Write position and value pairs
            for pos, count in zip(sorted_positions, sorted_counts):
                f.write(f"{pos}\t{count}\n")
        
        print(f"WIG file written successfully with {len(sorted_positions)} data points")
    else:
        print("No data to write to WIG file")


# Main loop: iterate through all regions as viewpoints
print(f"Starting virtual 4C analysis for {len(all_regions)} viewpoints...")
print(f"Resolution: {resolution}bp")

for viewpoint_name, viewpoint_start, viewpoint_end in all_regions:
    # Perform virtual 4C
    positions, counts, viewpoint_center, vp_start, vp_end = virtual_4c_for_viewpoint(
        viewpoint_name, viewpoint_start, viewpoint_end, 
        chr_name, resolution, query_start, query_end
    )
    
    # Create plot
    plot_virtual_4c(positions, counts, viewpoint_name, viewpoint_center, vp_start, vp_end,
                   chr_name, resolution, all_regions)
    
    # Write WIG file
    write_wig_file(positions, counts, viewpoint_name, vp_start, vp_end, chr_name, resolution)

print(f"\n{'='*80}")
print("All virtual 4C analyses complete!")
print(f"{'='*80}")
