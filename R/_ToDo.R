# To Do

# run tests in _ggkmcif3_test.r file and try to fix bugs, document persistent errors

# find call to xfun and update as required

# fix forestplotUV and forestplotMV to use the new uvsum2 and msummary functions
# update documentation

# run Build > check and fix bugs


data("pembrolizumab")
ggkmcif3(data=pembrolizumab,c("os_time","os_status"),"cohort")
