DIR          := models/SoftsusyNMSSM
MODNAME      := SoftsusyNMSSM
WITH_$(MODNAME) := yes

ifeq ($(WITH_SoftsusyMSSM),yes)
LIBSoftsusyNMSSM_SRC  := \
		$(DIR)/nmssmUtils.cpp \
		$(DIR)/nmssmsoftpars.cpp \
		$(DIR)/nmssmsoftsusy.cpp \
		$(DIR)/nmssmsusy.cpp \
		$(DIR)/nmssm1loop.f

ifneq ($(findstring two_scale,$(SOLVERS)),)
LIBSoftsusyNMSSM_SRC  += \
		$(DIR)/SoftsusyNMSSM_two_scale.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_sugra_constraint.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_susy_scale_constraint.cpp
endif
endif

LIBSoftsusyNMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSoftsusyNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSoftsusyNMSSM_SRC)))

LIBSoftsusyNMSSM_DEP  := \
		$(LIBSoftsusyNMSSM_OBJ:.o=.d)

LIBSoftsusyNMSSM      := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj

all-$(MODNAME): $(LIBSoftsusyNMSSM)

clean-$(MODNAME)-dep:
		-rm -f $(LIBSoftsusyNMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSoftsusyNMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBSoftsusyNMSSM)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSoftsusyNMSSM_DEP) $(LIBSoftsusyNMSSM_OBJ): CPPFLAGS += $(EIGENFLAGS)

$(LIBSoftsusyNMSSM): $(LIBSoftsusyNMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

ALLDEP += $(LIBSoftsusyNMSSM_DEP)
ALLLIB += $(LIBSoftsusyNMSSM)
