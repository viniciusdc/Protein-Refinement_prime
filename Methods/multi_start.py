

def multi_start(data):
    """prints/save the output and statistics from the refinement process when multi-start enable"""
    # check the data type:
    if type(data) != dict:
        print(":: data type object not match with dict structure!")
        print(":: The process was interrupted")
        return exit()
