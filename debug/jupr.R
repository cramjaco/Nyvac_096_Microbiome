cranPackages = c(
# to connect to jupyter
"repr",
"IRdisplay",
"evaluate",
"crayon",
"pbdZMQ",
"devtools",
"uuid",
"digest")

new.packages <- cranPackages[!(cranPackages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)