################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Alignment.cpp \
../src/BamParser.cpp \
../src/Parser.cpp \
../src/Sniffles.cpp 

OBJS += \
./src/Alignment.o \
./src/BamParser.o \
./src/Parser.o \
./src/Sniffles.o 

CPP_DEPS += \
./src/Alignment.d \
./src/BamParser.d \
./src/Parser.d \
./src/Sniffles.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/Users/fsedlaze/Documents/workspace/Sniffles/lib/bamtools/include/ -I/Users/fsedlaze/Documents/workspace/Sniffles/lib/tclap-1.2.1/include -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


