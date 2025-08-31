source env/bin/activate #or change to your env path
git pull upstream main
sudo cpufreq-set -c 3 -g performance
#asv does recent commits first, so by letting it run for 55min, and the cron to 1hr
#it will always keep up with new commits, but also process the backlog
timeout 3300s asv run -j 4 --show-stderr --cpu-affinity 3 --skip-existing ALL
sudo cpufreq-set -c 3 -g powersave
asv publish
asv preview
