BUILD_DIR ?= build
CMAKE_FLAGS ?= -DCMAKE_BUILD_TYPE=Release

.PHONY: all clean test install configure

all: configure
	cmake --build $(BUILD_DIR)

configure:
	cmake -B $(BUILD_DIR) $(CMAKE_FLAGS)

test: all
	cd $(BUILD_DIR) && ctest --output-on-failure

install: all
	cmake --install $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
