DIR_INC = ./include
DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_LIB = ./lib

CC = gcc
CFLAGS = -g -Wall -I $(DIR_INC)
LFLAGS = -lm

SRC = $(wildcard $(DIR_SRC)/*.c)
OBJ = $(patsubst %.c, $(DIR_OBJ)/%.o, $(notdir $(SRC)))

TAR = main

$(TAR):$(OBJ)
	$(CC) $^ -o $@ $(LFLAGS)
$(DIR_OBJ)/%.o:$(DIR_SRC)/%.c
	$(CC) $(CFLAGS) -c -fPIC $< -o $@
	
.PHONY:clean
clean:
	rm -rf $(OBJ)
	#rm -rf $(LIB_TAR)
	rm -rf $(TAR)

