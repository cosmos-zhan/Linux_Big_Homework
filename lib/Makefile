DIR_INC = ../include
DIR_SRC = ../src
DIR_OBJ = ../obj
DIR_LIB = ../lib

CC = gcc
CFLAGS = -g -Wall -I $(DIR_INC)
LFLAGS = -lm

SRC = $(wildcard $(DIR_SRC)/*.c)
OBJ = $(patsubst %.c, $(DIR_OBJ)/%.o, $(notdir $(SRC)))
LIB_OBJ = $(filter-out $(DIR_OBJ)/main.o, $(OBJ))
LIB_TAR = $(DIR_LIB)/libMethod.so

$(LIB_TAR):$(LIB_OBJ)
	$(CC) -shared $^ -o $@
$(DIR_OBJ)/%.o:$(DIR_SRC)/%.c
	$(CC) $(CFLAGS) -c -fPIC $< -o $@	
.PHONY:clean
clean:
	rm -rf $(OBJ)
	rm -rf $(LIB_TAR)
	#rm -rf $(TAR)
