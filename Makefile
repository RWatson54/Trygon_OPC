APP_NAME := trygon

OP2_DIR := /home/raw54/Source/OP2-DSL

APP_ENTRY := $(APP_NAME).cpp
APP_ENTRY_MPI := $(APP_ENTRY)

OP2_LIBS_WITH_HDF5 := true

VARIANT_FILTER_OUT := %vec

include $(OP2_DIR)/makefiles/common.mk
include $(OP2_DIR)/makefiles/c_app.mk

