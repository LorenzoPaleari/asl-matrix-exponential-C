import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.font_manager import FontProperties
import re

path = "../data_collection/"
versions = ["base_no_blas", "opt1_basic_no_blas", "opt3_basic_no_blas", "opt4_inline_basic_no_blas", "opt5_blocking_basic_no_blas", "opt6_vect_basic_no_blas"]
flops = []
runtime = []
performance = []
input_size = [4, 8, 16, 32, 64, 128, 256, 512, 1024]

for version in versions:
    file = open(path + version + "/flops_clean.txt", "r")
    lines = file.readlines()
    file.close()

    return_index = 6
    dimension_index = 0
    for line in lines[1:]:
        if len(line) < 15:
            dimension_index = 0
        numbers = re.findall(r'\b(\d+(?:\.\d+)?)N\b', line[15:])
        try:
            numbers.append((re.findall(r'\d+(?:\.\d+)?', line[15:])[-1:][0]))
            flops.append(float(numbers[0]) * (input_size[dimension_index]**3) + float(numbers[1]) * (input_size[dimension_index]**2) + float(numbers[2]) * input_size[dimension_index] + float(numbers[3]))
            dimension_index += 1
        except:
            pass

    for return_index in range(0, 6):
        for inputs in input_size:
            df = pd.read_csv(path + version + "/data/" + str(inputs) + ".csv")
            runtime.append(df[" Cycles"][return_index])

for i in range(0, len(runtime)):
    performance.append(flops[i] / runtime[i])


for return_index in range(0, 6):
# Set dimension
    fig, ax = plt.subplots(figsize=(14, 9))

    plt.rc('font', family='sans serif', size=12)
    font = FontProperties()
    font.set_family('sans serif')

    # Set Title
    title = ax.set_title('Intel(R) Core(TM) i5-6400 @ 2700Mhz \nL1:64KB L2:256KB L3:6MB \nCompiler: gcc 10.2.1', fontproperties=font, fontsize=16, fontweight='bold', pad=35, ha='left')
    title.set_position((-0.03, 1.4))


    # Create the line plot
    plt.xscale('log', base=2)

    ax.plot(input_size, performance[0 + return_index*9: 9 + return_index*9], color='black', linewidth=1.5, marker="o", markersize=10, label='Base no_blas',
        markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size, performance[54 + return_index*9: 63 + return_index*9], color='green', linewidth=1.5, marker="o", markersize=10, label='Opt1 no_blas',
        markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size, performance[108 + return_index*9: 117 + return_index*9], color='orange', linewidth=1.5, marker="o", markersize=10, label='Opt2 no_blas',
        markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size, performance[162 + return_index*9: 171 + return_index*9], color='blue', linewidth=1.5, marker="o", markersize=10, label='Inline no_blas',
        markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size, performance[216 + return_index*9: 225 + return_index*9], color='pink', linewidth=1.5, marker="o", markersize=10, label='Blocking no_blas',
        markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size, performance[270 + return_index*9: 279 + return_index*9], color='red', linewidth=3, marker="o", markersize=10, label='Vect no_blas',
        markeredgecolor='white', zorder=20, clip_on=False)



    # Set x axis label, limits and ticks
    ax.set_xlabel('Input Size', fontname='sans serif', fontsize=14, labelpad=4)
    ax.set_xlim(3, 1100)  # y limits from 0 to 1
    ticks = [4, 8, 16, 32, 64, 128, 256, 512, 1024]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["4", "8", "16", "32", "64", "128", "256", "512", "1024"])


    # Set y axis label, limits and ticks
    ax.set_ylabel('Performance [Flops/Cycle]', fontname='sans serif', fontsize=14, rotation='horizontal')
    ax.set_ylim(0, 6)  # y limits from 0 to 1
    ax.yaxis.grid(True, linewidth=2, c="white")
    ax.yaxis.set_label_coords(0.07, 1.02)

    # Tick and border settings
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params(axis='y', length=0, pad=10, labelsize=14)
    ax.tick_params(axis='x', pad=6, width=2, length=5, labelsize=14)

    #LEGEND
    #box0 = ax.get_position()
    #ax.set_position([box0.x0, box0.y0, box0.width * 0.94, box0.height])
    ax.legend(loc='upper right', bbox_to_anchor=(1, 1), facecolor='white')


    # Set the background color to grey
    ax.set_facecolor('#e5e5e5')
    fig.subplots_adjust(left=0.05, right=0.98, top=0.85, bottom=0.07)

    # Show the plot
    plt.show()
    fig.savefig('Images/base_No_blas' + str(6 - return_index) + '.svg', format='svg')


