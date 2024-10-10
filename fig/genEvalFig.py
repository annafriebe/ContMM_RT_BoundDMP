import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import csv
import json



def draw_bounds_8(axes, ax_row_ind, filename):
	dmp_bound_state_3 = []
	dmp_bound_overall = []
	x = []
	# the lines contain the accumulation index + 8 per-state betas + 
	# 8 per-state lower p_{wd} + 8 per-state higher p_{wd} + 8 per-state dmp bound +
	# the overall dmp bound
	dmp3_ind = 1 + 8*3 + 2
	dmp_overall_ind = 1 + 8*4
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter=',')
		for row in csv_reader:
			if len(row) > 2: 
				x.append(int(row[0]))
				dmp_bound_state_3.append(float(row[dmp3_ind]))
				dmp_bound_overall.append(float(row[dmp_overall_ind]))
	axes[ax_row_ind][0].plot(x, dmp_bound_state_3, color='black')
	axes[ax_row_ind][1].plot(x, dmp_bound_overall, color='black')

def draw_bounds_2(axes, ax_row_ind, filename):
	dmp_bound_worst = []
	dmp_bound_overall = []
	x = []
	# the lines contain the accumulation index + 2 per-state betas + 
	# 2 per-state lower p_{wd} + 2 per-state higher p_{wd} + 2 per-state dmp bound +
	# the overall dmp bound
	dmpworst_ind = 1 + 2*3 
	dmp_overall_ind = 1 + 2*4
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter=',')
		for row in csv_reader:
			if len(row) > 2: 
				x.append(int(row[0]))
				dmp_bound_worst.append(float(row[dmpworst_ind]))
				dmp_bound_overall.append(float(row[dmp_overall_ind]))
	axes[ax_row_ind][0].plot(x, dmp_bound_worst, color='green')
	axes[ax_row_ind][1].plot(x, dmp_bound_overall, color='green')

def index_from_Q_k(Q, k):
	if Q == 60000:
		if k == 10:
			return 0
		elif k == 8:
			return 1
	if Q == 70000:
		if k == 8:
			return 2
		elif k == 6:
			return 3
	if Q == 80000:
		if k == 8:
			return 4
		elif k == 6:
			return 5
	raise RuntimeError("index_from_Q_k: Unknown Q and k")


def draw_sim_cont(axes, filename):
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter=',')
		for row in csv_reader:
			if row[0].isnumeric():
				dmp_3 = float(row[4])
				dmp_all = float(row[5])
				index = index_from_Q_k(int(row[1]), int(row[3]))
				axes[index][0].hlines(y=dmp_3, xmin=1, xmax=10, color='red')
				axes[index][1].hlines(y=dmp_all, xmin=1, xmax=10, color='red')

def dm_sum_from_log(filename):
	with open(filename) as json_file:
		missed_deadlines = json.load(json_file)['missed deadlines']
		sum_missed_deadlines = 0
		for elem in missed_deadlines:
			if elem[0]>= 2000:
				sum_missed_deadlines += elem[1]
		return sum_missed_deadlines
		


def dmr_from_logs(Q, n, k):
	sum_missed_deadlines = 0
	for i in range(1, 11):
		filename = "../data/linux_dl_results/dl_" + str(Q) + "_" + str(n) + \
			"_" + str(k) + "_" + str(i) + ".json"
		sum_missed_deadlines += dm_sum_from_log(filename)
	return sum_missed_deadlines/480000

def dmp_from_prosit_log(filename, deadline):
	with open(filename) as prosit_file:
		for line in prosit_file:
			split_line = line.split()
			if split_line[0].isnumeric() and (int(split_line[0]) == deadline):
				return 1 - float(split_line[1])
	raise RuntimeError("dmp_from_prosit_log, deadline not found")

def dmp_ind(Q, k):
	filename = "../data/ind/result.csv"
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter=',')
		for row in csv_reader:
			if row[0].isnumeric():
				if (int(row[1]) == Q*1000) and (int(row[3]) == k):
					return float(row[4])
	raise RuntimeError("dmp_ind, Unknown Q and k")


fig, axs = plt.subplots(7, 2)
for i in range(6):
	axs[i][0].axis([1, 10, 0, 0.6])
	axs[i][1].axis([1, 10, 0, 0.01])

# Workaround, room for legend
axs[6, 0].axis('off')
axs[6, 1].axis('off')

dmp_bound_filenames = [\
"../data/dmp_bound/pzwl_dmp_control_bound_60000_5_10.csv",\
"../data/dmp_bound/pzwl_dmp_control_bound_60000_5_8.csv",\
"../data/dmp_bound/pzwl_dmp_control_bound_70000_4_8.csv",\
"../data/dmp_bound/pzwl_dmp_control_bound_70000_4_6.csv",\
"../data/dmp_bound/pzwl_dmp_control_bound_80000_4_8.csv",\
"../data/dmp_bound/pzwl_dmp_control_bound_80000_4_6.csv"\
]
for i in range(6):
	draw_bounds_8(axs, i, dmp_bound_filenames[i])

dmp_merged_bound_filenames = [\
"../data/dmp_bound/pzwl_dmp_control_merged_bound_60000_5_10.csv",\
"../data/dmp_bound/pzwl_dmp_control_merged_bound_60000_5_8.csv",\
"../data/dmp_bound/pzwl_dmp_control_merged_bound_70000_4_8.csv",\
"../data/dmp_bound/pzwl_dmp_control_merged_bound_70000_4_6.csv",\
"../data/dmp_bound/pzwl_dmp_control_merged_bound_80000_4_8.csv",\
"../data/dmp_bound/pzwl_dmp_control_merged_bound_80000_4_6.csv"\
]
for i in range(6):
	draw_bounds_2(axs, i, dmp_merged_bound_filenames[i])

sim_cont_filename = "../data/sim_cont/dmp_result.csv"	
draw_sim_cont(axs, sim_cont_filename)

Qs = [60, 70, 80]
ns = [5, 4, 4]
ks = [[10, 8], [8, 6], [8, 6]]

for i in range(3):
	for j in range(2):
		dmr_fig = dmr_from_logs(Qs[i], ns[i], ks[i][j])
		fig_ind = i*2 + j
		axs[fig_ind][1].hlines(y=dmr_fig, xmin=1, xmax=10, color='yellow')


for i in range(3):
	prosit_filename = "../data/prosit/result_" + str(Qs[i]) + "000_" + str(ns[i]) + ".txt"
	for j in range(2):
		deadline = ks[i][j]*200//(ns[i])
		dmp_prosit = dmp_from_prosit_log(prosit_filename, deadline)
		fig_ind = i*2 + j
		axs[fig_ind][1].hlines(y=dmp_prosit, xmin=1, xmax=10, color='black', linestyles='dashed')
		
for i in range(3):
	for j in range(2):
		dmp_fig = dmp_ind(Qs[i], ks[i][j])
		fig_ind = i*2 + j
		axs[fig_ind][1].hlines(y=dmp_fig, xmin=1, xmax=10, color='blue')
		
# titles
axs[0][0].set_title("Worst state")
axs[0][1].set_title("Overall")


axs[0][0].set_yticks(ticks = [0, 0.2, 0.4, 0.6], labels = ['0', '0.2', '0.4', ''])
axs[0][1].set_yticks(ticks = [0.002, 0.004, 0.006, 0.008, 0.01], \
	labels = ['', '', '', '', '0.01'])
axs[1][0].set_yticks(ticks = [0, 0.2, 0.4, 0.6], labels = ['0', '0.2', '0.4', ''])
axs[1][1].set_yticks(ticks = [0.002, 0.004, 0.006, 0.008, 0.01], \
	labels = ['', '', '', '', '0.01'])
axs[2][0].set_yticks(ticks = [0, 0.2, 0.4, 0.6], labels = ['0', '0.2', '0.4', ''])
axs[2][1].set_yticks(ticks = [0.002, 0.004, 0.006, 0.008, 0.01], \
	labels = ['', '', '', '', '0.01'])
axs[3][0].set_yticks(ticks = [0, 0.2, 0.4, 0.6], labels = ['0', '0.2', '0.4', ''])
axs[3][1].set_yticks(ticks = [0.002, 0.004, 0.006, 0.008, 0.01], \
	labels = ['', '', '', '', '0.01'])
axs[4][0].set_yticks(ticks = [0, 0.2, 0.4, 0.6], labels = ['0', '0.2', '0.4', ''])
axs[4][1].set_yticks(ticks = [0.002, 0.004, 0.006, 0.008, 0.01], \
	labels = ['', '', '', '', '0.01'])
axs[5][0].set_yticks(ticks = [0, 0.2, 0.4, 0.6], labels = ['0', '0.2', '0.4', ''])
axs[5][1].set_yticks(ticks = [0.002, 0.004, 0.006, 0.008, 0.01], \
	labels = ['', '', '', '', '0.01'])

axs[0][0].set_ylabel(".06/5/10 \n $p_{dm}$")
axs[1][0].set_ylabel(".06/5/8 \n $p_{dm}$")
axs[2][0].set_ylabel(".07/4/8 \n $p_{dm}$")
axs[3][0].set_ylabel(".07/4/6 \n $p_{dm}$")
axs[4][0].set_ylabel(".08/4/8 \n $p_{dm}$")
axs[5][0].set_ylabel(".08/4/6 \n $p_{dm}$")

axs[0][0].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['', '', '', '', ''])
axs[0][1].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['', '', '', '', ''])
axs[1][0].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['', '', '', '', ''])
axs[1][1].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['', '', '', '', ''])
axs[2][0].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['', '', '', '', ''])
axs[2][1].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['', '', '', '', ''])
axs[3][0].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['', '', '', '', ''])
axs[3][1].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['', '', '', '', ''])
axs[4][0].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['', '', '', '', ''])
axs[4][1].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['', '', '', '', ''])
axs[5][0].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['2', '4', '6', '8', '10'])
axs[5][1].set_xticks(ticks = [2, 4, 6, 8, 10], labels = ['2', '4', '6', '8', '10'])
axs[5][0].set_xlabel("Accumulation period")
axs[5][1].set_xlabel("Accumulation period")

axs[0][0].grid(visible=True)
axs[0][1].grid(visible=True)
axs[1][0].grid(visible=True)
axs[1][1].grid(visible=True)
axs[2][0].grid(visible=True)
axs[2][1].grid(visible=True)
axs[3][0].grid(visible=True)
axs[3][1].grid(visible=True)
axs[4][0].grid(visible=True)
axs[4][1].grid(visible=True)
axs[5][0].grid(visible=True)
axs[5][1].grid(visible=True)


# legends
#red_plot = mpatches.Patch(color="red", label = "Experimental")
black_plot = mlines.Line2D([],[],color="black", label = "Bound-8 on $p_{dm}$")
green_plot = mlines.Line2D([],[],color="green", label = "Bound-2 on $p_{dm}$")
yellow_plot = mlines.Line2D([],[],color="yellow", label = "Linux CBS experimental DMR")
dashed_plot = mlines.Line2D([],[],color="black", linestyle="dashed", label = "PROSIT-derived $p_{dm}$")
red_plot = mlines.Line2D([],[],color="red", label = "Sim-Cont $p_{dm}$ from simulation")
blue_plot = mlines.Line2D([],[],color="blue", label = "Ind $p_{dm}$ assuming independence")
#axs[5][0].legend(handles=[black_plot, yellow_plot, dashed_plot, red_plot, blue_plot], loc='upper left', )
#fig.legend(handles=[black_plot, green_plot, yellow_plot, dashed_plot, red_plot, blue_plot], \
#	loc = 'lower center', ncol = 2, fontsize='x-small')
fig.legend(handles=[black_plot, red_plot, blue_plot, green_plot, yellow_plot, dashed_plot], \
	loc = 'lower center', ncol = 2, fontsize='x-small')
#fig.legend(handles=[black_plot, yellow_plot, dashed_plot, red_plot, blue_plot], loc = 10)
fig.savefig("eval_fig.png")
fig.savefig("eval_fig.eps")


