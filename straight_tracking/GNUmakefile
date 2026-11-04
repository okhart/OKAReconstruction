
########################################################################
#                                                                      #
# Default target.                                                      #
#                                                                      #
########################################################################

.PHONY: default
default: tracking

.PHONY: clean
clean:
	rm -f $(exe_files) $(obj_files)

exe_files := tracking

obj_files := *.o #<-this-is-for-normal-RK-#


tracking: 
	g++ -I include *.cpp -o tracking `root-config --cflags --glibs` -fPIC



