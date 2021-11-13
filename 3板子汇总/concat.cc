#include <fstream>
#include <filesystem>
#include <string>

using namespace std;

int main() {
    ofstream fout("collection.md");
    for (auto subdir : filesystem::directory_iterator("template")) {
        fout << "\n# `" << subdir.path().filename() << "`\n\n";
        for (auto filename : filesystem::directory_iterator(subdir)) {
            if (filename.path().extension() == ".md") {
                fout << "\n## `" << filename.path().filename() << "`\n\n";
                ifstream fin(filename.path());
                string line;
                while (getline(fin, line)) fout << line << "\n";
                fout << "\n";
            }
            if (filename.path().extension() == ".cpp") {
                fout << "\n## `" << filename.path().filename() << "`\n\n";
                fout << "\n```cpp\n";
                ifstream fin(filename.path());
                string line;
                while (getline(fin, line)) fout << line << "\n";
                fout << "```\n\n";
            }
        }
    }
}
