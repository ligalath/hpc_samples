#define LOGD(content, ...) printf("==== djc debug ==== file: %s, line: %d, info: " content "\n",__FILE__, __LINE__, ##__VA_ARGS__)
#define LOGE(content, ...) printf("==== djc error ==== file: %s, line: %d, info: " content "\n",__FILE__, __LINE__, ##__VA_ARGS__)