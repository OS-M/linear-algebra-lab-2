import matplotlib.pyplot as plt
import numpy as np
import sys

lines = open(sys.argv[1], "r").read()
step = int(sys.argv[2])

plot_count = len(lines.strip().split("===")) - 1
fig = plt.figure(figsize=(5 * plot_count, 5))
plot_num = 1

for plot in lines.strip().split("==="):
    if len(plot) == 0:
        break
    plot = plot.strip().split('\n')
    name = plot[0]
    x_label = plot[1]
    y_label = plot[2]
    xs = [int(i) for i in plot[3].strip().split(' ')]
    xs = xs[::step]

    ax = fig.add_subplot(1, plot_count, plot_num)
    # ax.set_xticklabels(xs)
    ax.set_xticks(xs)
    ax.title.set_text(name)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    for i in range(4, len(plot), 2):
        line_name = plot[i]
        y = np.array([float(j) for j in plot[i + 1].strip().split(' ')])
        ys = []
        for j in range(0, len(y), step):
            ys.append(np.sum(y[j:j + step]))
        ax.bar(xs, ys, label=line_name, width=30)

    plot_num += 1
    plt.legend()

# plt.show()
plt.savefig(sys.argv[1][:-4] + ".png")
