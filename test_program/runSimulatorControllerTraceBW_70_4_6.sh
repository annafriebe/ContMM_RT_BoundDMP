# !/bin/bash

echo  performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

# disable admission control for core pinning (for CBS)
echo -1 > /proc/sys/kernel/sched_rt_runtime_us
mkdir /sys/fs/cgroup/cpuset/Set2
echo 2 > /sys/fs/cgroup/cpuset/Set2/cpuset.cpus
echo 0 > /sys/fs/cgroup/cpuset/Set2/cpuset.mems
echo 0 > /sys/fs/cgroup/cpuset/Set2/cpuset.sched_load_balance
echo 1 > /sys/fs/cgroup/cpuset/Set2/cpuset.cpu_exclusive

mkdir /sys/fs/cgroup/cpuset/Set1
echo 1 > /sys/fs/cgroup/cpuset/Set1/cpuset.cpus
echo 0 > /sys/fs/cgroup/cpuset/Set1/cpuset.mems
echo 0 > /sys/fs/cgroup/cpuset/Set1/cpuset.sched_load_balance
echo 1 > /sys/fs/cgroup/cpuset/Set1/cpuset.cpu_exclusive


# disable USB
echo '1-1' | tee /sys/bus/usb/drivers/usb/unbind

./simulator/ripRTSimulate 1 &
trace-cmd record -e sched -o ../data/traces_reports_csv/controlTraceDeadlineBW.dat ./controller/ripControlDeadlineBW 70000 500000 4 6 "../data/linux_dl_results/dl_70_4_6_1.json" 1

echo  ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

# cool down, ensure simulation stopped
sleep 30s

echo  performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

./simulator/ripRTSimulate 1 &
trace-cmd record -e sched -o ../data/traces_reports_csv/controlTraceDeadlineBW.dat ./controller/ripControlDeadlineBW 70000 500000 4 6 "../data/linux_dl_results/dl_70_4_6_2.json" 1

echo  ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor
# cool down, ensure simulation stopped
sleep 30s

echo  performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

./simulator/ripRTSimulate 1 &
trace-cmd record -e sched -o ../data/traces_reports_csv/controlTraceDeadlineBW.dat ./controller/ripControlDeadlineBW 70000 500000 4 6 "../data/linux_dl_results/dl_70_4_6_3.json" 1

echo  ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor
# cool down, ensure simulation stopped
sleep 30s

echo  performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

./simulator/ripRTSimulate 1 &
trace-cmd record -e sched -o ../data/traces_reports_csv/controlTraceDeadlineBW.dat ./controller/ripControlDeadlineBW 70000 500000 4 6 "../data/linux_dl_results/dl_70_4_6_4.json" 1

echo  ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor
# cool down, ensure simulation stopped
sleep 30s

echo  performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

./simulator/ripRTSimulate 1 &
trace-cmd record -e sched -o ../data/traces_reports_csv/controlTraceDeadlineBW.dat ./controller/ripControlDeadlineBW 70000 500000 4 6 "../data/linux_dl_results/dl_70_4_6_5.json" 1

echo  ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor
# cool down, ensure simulation stopped
sleep 30s

echo  performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

./simulator/ripRTSimulate 1 &
trace-cmd record -e sched -o ../data/traces_reports_csv/controlTraceDeadlineBW.dat ./controller/ripControlDeadlineBW 70000 500000 4 6 "../data/linux_dl_results/dl_70_4_6_6.json" 1

echo  ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor
# cool down, ensure simulation stopped
sleep 30s

echo  performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

./simulator/ripRTSimulate 1 &
trace-cmd record -e sched -o ../data/traces_reports_csv/controlTraceDeadlineBW.dat ./controller/ripControlDeadlineBW 70000 500000 4 6 "../data/linux_dl_results/dl_70_4_6_7.json" 1

echo  ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor
# cool down, ensure simulation stopped
sleep 30s

echo  performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

./simulator/ripRTSimulate 1 &
trace-cmd record -e sched -o ../data/traces_reports_csv/controlTraceDeadlineBW.dat ./controller/ripControlDeadlineBW 70000 500000 4 6 "../data/linux_dl_results/dl_70_4_6_8.json" 1

echo  ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor
# cool down, ensure simulation stopped
sleep 30s

echo  performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

./simulator/ripRTSimulate 1 &
trace-cmd record -e sched -o ../data/traces_reports_csv/controlTraceDeadlineBW.dat ./controller/ripControlDeadlineBW 70000 500000 4 6 "../data/linux_dl_results/dl_70_4_6_9.json" 1

echo  ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor
# cool down, ensure simulation stopped
sleep 30s

echo  performance > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  performance > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

./simulator/ripRTSimulate 1 &
trace-cmd record -e sched -o ../data/traces_reports_csv/controlTraceDeadlineBW.dat ./controller/ripControlDeadlineBW 70000 500000 4 6 "../data/linux_dl_results/dl_70_4_6_10.json" 1

echo  ondemand > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu1/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor
echo  ondemand > /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor

#reenable USB
echo '1-1' | tee /sys/bus/usb/drivers/usb/bind
