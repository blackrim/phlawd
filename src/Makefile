CC = g++
CFLAGS = -msse3 -DSSE2 -O3 -g -w -Wall -fopenmp

CPP_SRCS += \
sequence.cpp \
tree.cpp \
tree_reader.cpp \
node.cpp \
tree_utils.cpp \
fasta_util.cpp \
DBSeq.cpp \
SQLiteConstructor.cpp \
SQLiteProfiler.cpp \
SQLiteTreeNameConvertor.cpp \
SQLiteDBController.cpp \
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
sequence.o \
tree.o \
tree_reader.o \
node.o \
tree_utils.o \
fasta_util.o \
DBSeq.o \
SQLiteConstructor.o \
SQLiteProfiler.o \
SQLiteTreeNameConvertor.o \
SQLiteDBController.o \
SWPS3_DynProgr_scalar.o \
SWPS3_DynProgr_sse_byte.o \
SWPS3_DynProgr_sse_short.o \
SWPS3_Same_seq_pthread.o \
SWPS3_debug.o \
SWPS3_fasta.o \
SWPS3_matrix.o \
SWPS3_swps3.o \
Same_seq_pthread.o \
SmithWatermanGotoh.o \
main.o \
utils.o 

CPP_DEPS += \
sequence.d \
tree.d \
tree_reader.d \
node.d \
tree_utils.d \
fasta_util.d \
DBSeq.d \
SQLiteConstructor.d \
SQLiteProfiler.d \
SQLiteTreeNameConvertor.d \
SQLiteDBController.d \
SWPS3_DynProgr_scalar.d \
SWPS3_DynProgr_sse_byte.d \
SWPS3_DynProgr_sse_short.d \
SWPS3_Same_seq_pthread.d \
SWPS3_debug.d \
SWPS3_fasta.d \
SWPS3_matrix.d \
SWPS3_swps3.d \
Same_seq_pthread.d \
SmithWatermanGotoh.d \
main.d \
utils.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.cpp
	$(CC) $(CFLAGS)  -L../deps/linux/ -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"

LIBS := -lm -lsqlitewrapped -lsqlite3 -fopenmp

RM := rm -rf

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: PHLAWD

# Tool invocations
PHLAWD: $(OBJS)
	@echo 'Building target: $@'
	$(CC) -msse3 -DSSE2 -L../deps/linux -o "PHLAWD" $(OBJS) $(LIBS) -lpthread
	@echo 'Finished building target: $@'
	@echo ' '

clean:
	-$(RM) *.o *.d
	-@echo ' '
