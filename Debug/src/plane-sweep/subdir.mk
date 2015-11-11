################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/plane-sweep/Main.cpp \
../src/plane-sweep/MyHeap.cpp \
../src/plane-sweep/MyList.cpp \
../src/plane-sweep/PlaneSweep.cpp 

OBJS += \
./src/plane-sweep/Main.o \
./src/plane-sweep/MyHeap.o \
./src/plane-sweep/MyList.o \
./src/plane-sweep/PlaneSweep.o 

CPP_DEPS += \
./src/plane-sweep/Main.d \
./src/plane-sweep/MyHeap.d \
./src/plane-sweep/MyList.d \
./src/plane-sweep/PlaneSweep.d 


# Each subdirectory must supply rules for building sources it contributes
src/plane-sweep/%.o: ../src/plane-sweep/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/Users/fsedlaze/Documents/workspace/Sniffles/lib/bamtools/include/ -I/Users/fsedlaze/Documents/workspace/Sniffles/lib/tclap-1.2.1/include -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


