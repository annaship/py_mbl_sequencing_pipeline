import logging
import logging.handlers
import os

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

# Remove the default FileHandlers if present. 
handler = logging.handlers.RotatingFileHandler(os.path.join(os.getcwd(),"pipeline.log"), maxBytes=10000000, backupCount=10)
logger.addHandler(handler)
handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")) 


