import logging
import logging.handlers
import os
import sys
import collections

# create logger with 'spam_application'
logger = logging.getLogger('')
# create console handler with a higher log level
#ch = logging.StreamHandler()
# create formatter and add it to the handlers
#formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
#ch.setFormatter(formatter)
# add the handlers to the logger
#logger.addHandler(ch)

# this gets rid of double messages:
logging.basicConfig(level=logging.DEBUG,format=FORMAT)

running_command_line = " ".join(sys.argv)
arg_names = sys.argv[1::2]
arg_vals = sys.argv[2::2]
args = dict(zip(arg_names, arg_vals))

logger.info("START\n%s", running_command_line)
# logger.debug("Number of arguments: ", len(sys.argv))
# logger.info("The arguments are: " , str(sys.argv))

log_file_name = "pipeline_%s_%s_%s_%s.log" % (args['-r'], args['-p'], args['-lane_name'], args['-db_name'])
# Remove the default FileHandlers if present.
handler = logging.handlers.RotatingFileHandler(os.path.join(os.getcwd(), log_file_name), maxBytes=10000000, backupCount=10)
logger.addHandler(handler)
logger.debug("log path: %s" % handler.baseFilename)
handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))


