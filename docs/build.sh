#/bin/bash

# Jupyter-build doesn't have an option to automatically show the 
# saved reports, which makes it difficult to debug the reasons for 
# build failures in CI. This is a simple wrapper to handle that.

REPORTDIR=_build/html/reports

# TODO add -n for nitpick mode.
jupyter-book build -W --keep-going .
RETVAL=$?
if [ $RETVAL -ne 0 ]; then
    if [ -e $REPORTDIR ]; then
      echo "Error occured; showing saved reports"
      cat $REPORTDIR/*
    fi
else
    # Clear out any old reports
    rm -f $REPORTDIR/*
fi
exit $RETVAL
