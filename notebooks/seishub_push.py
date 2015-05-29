

#! /usr/bin/env python

import os
from obspy.xseed import Parser

#very messisly created string containing the seishub push command:

check = False  
while check == False:
    
    callfile = raw_input("enter file name to push to SeisHub, e.g. event1.xml:  ")
    
    curl = 'curl -v --data-binary '
    prefix = '@'
    file_name = prefix + callfile
    post = ' -u admin:admin -X POST '
    url_instance = 'http://localhost:8080/xml/seismology/event/'
    URL = url_instance + callfile
#will need to change http path to running instance of seishub
#combine components defined above to create command string
    seishub_push = curl + file_name + post + URL
  
    if ".xml" in callfile:       
        os.system(seishub_push)
#        print seishub_push
        check = True
        
    elif ".resp" in callfile:
        print "haven't yet incorperated RESP files into seishub push. Please use different format "
        check = False  
        
    elif "seed" in callfile:
#won't work if the word 'seed' is in the filename
        check2 = False
        while check2 == False:

            convertfile = raw_input("do you wish to convert to XML before pushing?  ")    
            if 'yes' in convertfile:
            #conversion from XMl to mSEED - not working yet
#                sp = Parser(file_name)
#                sp.writeXSEED(file_name)
                print 'converted to .xml'
                print seishub_push
#                os.system(seishub_push)        
                check2 = True
            #need to use new .XML name here for the push, as oppossed to originally inputted name       

            elif 'no' in convertfile:  
                print seishub_push
                check2 = True
        del check2
        
    else:    
        print "response not recognised, please answer only yes/no "
        check = False

        
        check = True
del check

