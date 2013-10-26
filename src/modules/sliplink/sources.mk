include $(SRC_DIR)/modules/sliplink/mcMoves/sources.mk
include $(SRC_DIR)/modules/sliplink/diagnostics/sources.mk

modules_sliplink_= \
    $(modules_sliplink_mcMoves_) \
    $(modules_sliplink_diagnostics_) \
    modules/sliplink/SliplinkMcModule.cpp 

modules_sliplink_SRCS=$(addprefix $(SRC_DIR)/, $(modules_sliplink_))
modules_sliplink_OBJS=$(addprefix $(OBJ_DIR)/, $(modules_sliplink_:.cpp=.o))

