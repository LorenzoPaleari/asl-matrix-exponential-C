import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.font_manager import FontProperties
import re

path = "../data_collection/"
versions = ["base", "opt3_basic", "opt4_inline_basic", "opt5_blocking_basic", "opt6_vect_basic", "opt6_vect_lu", "lu"]
flops = []
runtime = []
performance = []
input_size = [4, 8, 16, 32, 64, 128, 256, 512, 1024]
scipy = [0.000000, 0.0, 0.105, 0.568, 0.0, 0.0, 0.110, 0.578, 0.0, 0.0, 0.086, 0.45, 0.0, 0.0, 0.033, 0.230,0.0,0.0, 0.024, 0.129, 0.0, 0.0, 0.020, 0.123]
for i in range(0, len(scipy)):
    scipy[i] *= 1000

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
            runtime.append(df[" Cycles"][return_index] / 2700000.0)

for i in range(0, len(runtime)):
    performance.append(flops[i] / runtime[i])

print(max(runtime))
print(runtime)

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

    ax.plot(input_size[5:], runtime[5 + return_index*9: 9 + return_index*9], color='black', linewidth=1.5, marker="o", markersize=10, label='Base',
        markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size[5:], runtime[59 + return_index*9: 63 + return_index*9], color='#a1a1a1', linewidth=1.5, marker="o", markersize=10, label='Base_optimized',
        markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size[5:], runtime[113 + return_index*9: 117 + return_index*9], color='#c3c3c3', linewidth=1.5, marker="o", markersize=10, label='Base_inlined',
        markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size[5:], runtime[167 + return_index*9: 171 + return_index*9], color='#b4b4b4', linewidth=1.5, marker="o", markersize=10, label='Base_Blocked',
         markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size[5:], runtime[221 + return_index*9: 225 + return_index*9], color='yellow', linewidth=1.5, marker="o", markersize=10, label='Base_vectorized',
         markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size[5:], runtime[275 + return_index*9: 279 + return_index*9], color='red', linewidth=3, marker="o", markersize=10, label='Lu_vectorized',
         markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size[5:], runtime[329 + return_index*9: 333 + return_index*9], color='green', linewidth=1.5, marker="o", markersize=10, label='Lu',
         markeredgecolor='white', zorder=20, clip_on=False)
    ax.plot(input_size[5:], scipy[0 + return_index*4: 4 + return_index*4], color='orange', linewidth=3, marker="o", markersize=10, label='Scipy',
         markeredgecolor='white', zorder=20, clip_on=False)



    # Set x axis label, limits and ticks
    ax.set_xlabel('Input Size', fontname='sans serif', fontsize=14, labelpad=4)
    ax.set_xlim(120, 1100)  # y limits from 0 to 1
    ticks = [128, 256, 512, 1024]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["128", "256", "512", "1024"])


    # Set y axis label, limits and ticks
    ax.set_ylabel('Runtime [Millisecond]', fontname='sans serif', fontsize=14, rotation='horizontal')
    ax.set_ylim(0, 2000)  # y limits from 0 to 1
    ax.yaxis.grid(True, linewidth=2, c="white")
    ax.yaxis.set_label_coords(0.052, 1.02)

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
    ax.legend(loc='upper left', bbox_to_anchor=(0, 1), facecolor='white')


    # Set the background color to grey
    ax.set_facecolor('#e5e5e5')
    fig.subplots_adjust(left=0.05, right=0.98, top=0.85, bottom=0.07)

    # Show the plot
    plt.show()
    fig.savefig('Images/runtime_' + str(6 - return_index) + '.svg', format='svg')


