DIR          := src
MODNAME      := libflexisusy

LIBFLEXI_SRC := \
		$(DIR)/coupling_monitor.cpp \
		$(DIR)/def.cpp \
		$(DIR)/dilog.f \
		$(DIR)/linalg.cpp \
		$(DIR)/lowe.cpp \
		$(DIR)/numerics.cpp \
		$(DIR)/rge.cpp \
		$(DIR)/stopwatch.cpp \
		$(DIR)/utils.cpp

ifneq ($(findstring two_scale,$(ALGORITMS)),)
LIBFLEXI_SRC += \
		$(DIR)/two_scale_composite_convergence_tester.cpp \
		$(DIR)/two_scale_convergence_tester.cpp \
		$(DIR)/two_scale_running_precision.cpp \
		$(DIR)/two_scale_solver.cpp
endif

LIBFLEXI_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFLEXI_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFLEXI_SRC)))

LIBFLEXI_DEP := \
		$(LIBFLEXI_OBJ:.o=.d)

LIBFLEXI     := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBFLEXI)

clean-$(MODNAME):
		rm -rf $(LIBFLEXI_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBFLEXI_DEP)
		rm -rf $(LIBFLEXI)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBFLEXI): $(LIBFLEXI_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBFLEXI_DEP)
ALLLIB += $(LIBFLEXI)
