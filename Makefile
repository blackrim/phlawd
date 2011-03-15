CC = g++
CFLAGS = -msse3 -DSSE2 -O3 -g -w -Wall

CPP_SRCS += \
DBSeq.cpp \
SQLiteConstructor.cpp \
SQLiteProfiler.cpp \
SQLiteTreeNameConvertor.cpp \
SWPS3_DynProgr_scalar.cpp \
SWPS3_DynProgr_sse_byte.cpp \
SWPS3_DynProgr_sse_short.cpp \
SWPS3_Same_seq_pthread.cpp \
SWPS3_debug.cpp \
SWPS3_fasta.cpp \
SWPS3_matrix.cpp \
SWPS3_swps3.cpp \
Same_seq_pthread.cpp \
SmithWatermanGotoh.cpp \
main.cpp \
utils.cpp 

OBJS += \
./DBSeq.o \
./SQLiteConstructor.o \
./SQLiteProfiler.o \
./SQLiteTreeNameConvertor.o \
./SWPS3_DynProgr_scalar.o \
./SWPS3_DynProgr_sse_byte.o \
./SWPS3_DynProgr_sse_short.o \
./SWPS3_Same_seq_pthread.o \
./SWPS3_debug.o \
./SWPS3_fasta.o \
./SWPS3_matrix.o \
./SWPS3_swps3.o \
./Same_seq_pthread.o \
./SmithWatermanGotoh.o \
./main.o \
./utils.o 

CPP_DEPS += \
./DBSeq.d \
./SQLiteConstructor.d \
./SQLiteProfiler.d \
./SQLiteTreeNameConvertor.d \
./SWPS3_DynProgr_scalar.d \
./SWPS3_DynProgr_sse_byte.d \
./SWPS3_DynProgr_sse_short.d \
./SWPS3_Same_seq_pthread.d \
./SWPS3_debug.d \
./SWPS3_fasta.d \
./SWPS3_matrix.d \
./SWPS3_swps3.d \
./Same_seq_pthread.d \
./SmithWatermanGotoh.d \
./main.d \
./utils.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: %.cpp
	$(CC) $(CFLAGS)  -Ldeps/linux/ -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"

LIBS := -lm -lsqlitewrapped -lsqlite3

RM := rm -rf

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: PHLAWD

# Tool invocations
PHLAWD: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: C++ Linker'
	$(CC) -msse3 -DSSE2 -Ldeps/linux -o "PHLAWD" $(OBJS) $(USER_OBJS) $(LIBS) -lpthread
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(EXECUTABLES)$(CC_DEPS)$(C_DEPS)$(CPP_DEPS)$(CXX_DEPS)$(C_UPPER_DEPS) PHLAWD
	-@echo ' '
