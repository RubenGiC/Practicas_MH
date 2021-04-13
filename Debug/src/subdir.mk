################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/PAR.cpp \
../src/Practicas_MH.cpp \
../src/random.cpp \
../src/result_algorithms.cpp 

OBJS += \
./src/PAR.o \
./src/Practicas_MH.o \
./src/random.o \
./src/result_algorithms.o 

CPP_DEPS += \
./src/PAR.d \
./src/Practicas_MH.d \
./src/random.d \
./src/result_algorithms.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


