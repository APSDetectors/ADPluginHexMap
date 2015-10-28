TOP=../..
include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE
#=============================

# The following gets rid of the -fno-implicit-templates flag on vxWorks, 
# so we get automatic template instantiation.
# This is what we want for the source files that define template functions.
ifeq (vxWorks,$(findstring vxWorks, $(T_A)))
CODE_CXXFLAGS=
endif

DBD += NDPluginHexMap.dbd

INC += NDPluginHexMap.h

LIBRARY_IOC += NDPluginHexMap

NDPluginHexMap_SRCS += NDPluginHexMap.cpp

include $(ADCORE)/ADApp/commonLibraryMakefile

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE

