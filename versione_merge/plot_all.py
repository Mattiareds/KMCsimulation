import numpy as np
import matplotlib.pyplot as plt
import os

# ============================================================
# 1) DATA FILES
# ============================================================
n_runs = 100
base_name = "run_"
occupation_file = "occupation.out"

folders = [f"{base_name}{i:03d}" for i in range(1, n_runs + 1)]

# ============================================================
# 2) NAMES TO WRITE — fill these in
# ============================================================

# --- Axis labels ---
x_label = "Total occupied sites"   # <-- write the x-axis label here
y_label = "Occupied sites"   # <-- write the y-axis label here
time_label = "Time of the simulation (s) "     # <-- write the top x-axis (time) label here

# --- Plot title (optional) ---
plot_title = ""   # <-- write a title here, or leave empty

# --- Legend names, one per column (in order y1..y6) ---
legend_names = [
    "1B",   # <-- name for column 1 (y1)
    "1C",   # <-- name for column 2 (y2)
    "0B-Pt",   # <-- name for column 3 (y3)
    "0C-Pt",   # <-- name for column 4 (y4)
    "0B-Au",   # <-- name for column 5 (y5)
    "0C-Au",   # <-- name for column 6 (y6)
]

# Maximum available value for columns 3, 4, 5, 6 (in this order)
max_values = [41, 64, 32, 49]   # <-- fill in the 4 max values here

# ============================================================
# 3) COLORS — one per column, change if you like
# ============================================================
colors = ["tab:blue", "tab:orange", "tab:green",
          "tab:red", "tab:purple", "tab:brown"]

# Show 1 tick out of every N on the top (time) axis,
# increase this if the time labels overlap.
time_tick_step = 10

# ============================================================
# 4) LOAD occupation.out FROM ALL RUNS, COMPUTE MEAN & STD
# ============================================================
all_y = []   # will become shape (n_valid_runs, n_points, 6)
x_ref = None

for folder in folders:
    path = os.path.join(folder, occupation_file)
    if not os.path.isfile(path):
        print(f"Warning: {path} not found, skipping.")
        continue

    data = np.loadtxt(path)
    x = data[:, 0]
    y = data[:, 1:7]

    if len(y) < 180:
        print(f"Warning: {path} has only {len(y)} rows (<175), excluding.")
        continue

    if x_ref is None:
        x_ref = x   # assume x values are the same across all runs

    all_y.append(y)

min_len = min(len(y) for y in all_y)
all_y = [y[:min_len] for y in all_y]
x_ref = x_ref[:min_len]

all_y = np.array(all_y)              # (n_valid_runs, n_points, 6)
mean_y = np.mean(all_y, axis=0)      # (n_points, 6)
std_y = np.std(all_y, axis=0)        # (n_points, 6)

x = x_ref

# ============================================================
# 5) LOAD TIME FROM occ####.out, AVERAGED OVER THE 100 RUNS
# ============================================================
# Assumption: occ####.out lives in the same run_* folder as
# occupation.out, and #### is the x value zero-padded to 4 digits.
time_marker = "time of the simulation:"

mean_time = []

for xi in x:
    suffix = f"{int(round(xi)):04d}"
    filename = f"occ{suffix}.out"

    times_for_this_x = []
    for folder in folders:
        path = os.path.join(folder, filename)
        if not os.path.isfile(path):
            continue
        with open(path) as f:
            lines = f.readlines()
        if len(lines) < 2:
            continue
        line = lines[1]   # second row of the file
        idx = line.find(time_marker)
        if idx == -1:
            continue
        value_str = line[idx + len(time_marker):].strip()
        try:
            times_for_this_x.append(float(value_str))
        except ValueError:
            continue

    if times_for_this_x:
        mean_time.append(np.mean(times_for_this_x))
    else:
        mean_time.append(np.nan)

mean_time = np.array(mean_time)

# ============================================================
# 6) PLOT
# ============================================================
fig, ax = plt.subplots(figsize=(8, 6))

for j, (color, name) in enumerate(zip(colors, legend_names)):
    ax.plot(x, mean_y[:, j], marker=".", linestyle="None", color=color, label=name)
    ax.fill_between(x, mean_y[:, j] - std_y[:, j], mean_y[:, j] + std_y[:, j],
                     color=color, alpha=0.2)

for j, m in zip(range(2, 6), max_values):
    ax.axhline(y=m, color=colors[j], linestyle="--", alpha=0.7)

ax.set_xlabel(x_label)
ax.set_ylabel(y_label)
ax.set_title(plot_title)
ax.legend()
ax.grid(True, alpha=0.3)

# --- top axis: time, aligned with the bottom x tick positions ---
ax_time = ax.twiny()
ax_time.set_xlim(ax.get_xlim())
ax_time.set_xticks(x[::time_tick_step])
ax_time.set_xticklabels([f"{t:.3g}" for t in mean_time[::time_tick_step]],
                         rotation=45, ha="left")
ax_time.set_xlabel(time_label)

plt.tight_layout()
plt.savefig("plot.png", dpi=300)
plt.show()
