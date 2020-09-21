source env/bin/activate
git pull upstream main
sudo cpufreq-set -c 3 -g performance
#asv does recent commits first, so by letting it run for 55min, and the cron to 1hr 
#it will always keep up with new commits, but also process the backlog
timeout 3300s asv run -j 4 --show-stderr --cpu-affinity 3 --skip-existing ALL
sudo cpufreq-set -c 3 -g powersave
asv publish
cd .asv
rm -rf docs
mv html docs #Github doesn't let you choose arbitary subfolders to serve so we have to use docs
rm -rf .git
git init
git remote add origin git@github.com:tskit-dev/msprime-asv.git
git checkout -b main
git add *
git commit -m "Automated asv upload"
git push -f origin main
cd ..
