# from matplotlib_venn import venn3
from matplotlib_venn import venn2
import argparse

import matplotlib.pyplot as plt

def plot_venn(AB, Ab, aB, output, plot_title, title_Ab, title_aB):
    plt.figure(figsize=(8, 8))

    print('AB:', AB)
    print('Ab:', Ab)
    print('aB:', aB)

    # Create scaled subsets for the venn diagram
    scaling_factor = 1000
    scaled_AB = AB / scaling_factor
    scaled_Ab = Ab / scaling_factor
    scaled_aB = aB / scaling_factor

    # Create a venn diagram scaled to the number of elements in each set
    # venn = venn2(subsets=(AB, Ab, aB), set_labels=(title_Ab, title_aB))
    venn = venn2(subsets=(scaled_Ab, scaled_aB, scaled_AB), set_labels=(title_Ab, title_aB))

    # Update the labels to reflect the actual counts
    venn.get_label_by_id('10').set_text(str(Ab))
    venn.get_label_by_id('01').set_text(str(aB))
    venn.get_label_by_id('11').set_text(str(AB))

    # Update the title
    # plt.title("contextsv and " + title_aB + " venn diagram (all SV types)")
    plt.title(plot_title)
    plt.savefig(output)
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a Venn diagram.')
    parser.add_argument('-a', type=int, required=True, help='Shared count')
    parser.add_argument('-b', type=int, required=True, help='False positive count')
    parser.add_argument('-c', type=int, required=True, help='False negative count')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file path')
    parser.add_argument('-a_title', type=str, required=True, help='Title for set A')
    parser.add_argument('-b_title', type=str, required=True, help='Title for set B')
    parser.add_argument('-c_title', type=str, required=True, help='Title for set C')

    args = parser.parse_args()

    plot_venn(args.a, args.b, args.c, args.output, args.a_title, args.b_title, args.c_title)
    print(f'Venn diagram saved to {args.output}')
