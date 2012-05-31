

class ConfigurationException(Exception):
    def __init__(self, value):
        super(ConfigurationException, self).__init__(value)
        self.value = value
    def __str__(self):
        return repr(self.value)    
