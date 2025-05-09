include make.inc
# Makefile
# 目标文件存放目录
BUILD_DIR= build

# 源文件目录
SRC_DIR := src

# 源文件列表
SRC := $(wildcard $(SRC_DIR)/*.cpp)

# 目标文件列表
OBJ := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC))

# 可执行文件名称
TARGET := main

# 默认目标
all: $(TARGET)

# 生成目标文件
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS) $(LIBS) $(INCLUDES)

# 生成可执行文件
$(TARGET): $(OBJ)
	$(CC) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS) $(INCLUDES)

# 清理生成的文件
clean:
	rm -rf $(BUILD_DIR) $(TARGET)