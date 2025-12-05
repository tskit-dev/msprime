source env/bin/activate #or change to your env path
git pull upstream main
sudo cpufreq-set -c 3 -g performance
#asv does recent commits first, so by letting it run for 55min, and the cron to 1hr
#it will always keep up with new commits, but also process the backlog
timeout 3300s asv run -j 4 --show-stderr --cpu-affinity 3 --skip-existing ALL

#instead you can test for a set of versions: uncomment the next line and comment line 6 to do so
#you can add more version to benshmark to the file version_hashes.txt
#asv run -j 4 --show-stderr --cpu-affinity 3 --skip-existing HASHFILE:benchmarks/version_hashes.txt
sudo cpufreq-set -c 3 -g powersave
asv publish
asv preview
