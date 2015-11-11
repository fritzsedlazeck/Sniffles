################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/tree/IntervallTree.cpp \
../src/tree/Scapegoat.cpp 

OBJS += \
./src/tree/IntervallTree.o \
./src/tree/Scapegoat.o 

CPP_DEPS += \
./src/tree/IntervallTree.d \
./src/tree/Scapegoat.d 


# Each subdirectory must supply rules for building sources it contributes
src/tree/%.o: ../src/tree/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/Users/fsedlaze/Documents/workspace/Sniffles/lib/bamtools/include/ -I/Users/fsedlaze/Documents/workspace/Sniffles/lib/tclap-1.2.1/include -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


