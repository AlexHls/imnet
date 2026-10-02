#include "file_io.h"
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <stdexcept>

void require(bool condition) {
  if (!condition) throw std::runtime_error("atomic save self-check failed");
}

int main() {
  namespace fs = std::filesystem;
  const auto directory = fs::temp_directory_path() /
      ("imnet-save-tests-" + std::to_string(std::random_device{}()));
  require(fs::create_directory(directory));
  const auto target = directory / "result.txt";
  const auto read = [&] {
    std::ifstream in(target);
    return std::string(std::istreambuf_iterator<char>(in), {});
  };
  using imyann::write_file_atomic;
  require(write_file_atomic(target.string(), [](auto &out) { out << "original"; }));
  const auto permissions = fs::status(target).permissions();
  for (bool stream_failure : {false, true}) {
    require(!write_file_atomic(target.string(), [&](auto &out) {
      out << "partial replacement";
      if (stream_failure) out.setstate(std::ios::badbit);
      else throw std::runtime_error("injected writer failure");
    }));
    require(read() == "original");
    require(std::distance(fs::directory_iterator(directory), fs::directory_iterator{}) == 1);
  }
  const auto link = directory / "linked.txt";
  fs::create_symlink(target, link);
  require(write_file_atomic(link.string(), [](auto &out) { out << "replacement"; }));
  require(fs::is_symlink(link) && read() == "replacement");
  require(fs::status(target).permissions() == permissions);
  require(!write_file_atomic(directory.string(), [](auto &out) { out << "bad"; }));
  require(!write_file_atomic((directory / "missing" / "file").string(), [](auto &) {}));
  // A destination becoming a directory before replacement must fail cleanly.
  const auto collision = directory / "collision";
  require(!write_file_atomic(collision.string(), [&](auto &out) {
    out << "complete";
    fs::create_directory(collision);
  }));
  require(fs::is_directory(collision) && read() == "replacement");
  require(std::distance(fs::directory_iterator(directory), fs::directory_iterator{}) == 3);
  require(imyann::save_trajectory(target.string(), {0, 1}, {1e4, 1e4}, {2e8, 2e8}));
  std::vector<double> times, rhos, temps;
  require(imyann::load_trajectory(target.string(), times, rhos, temps));
  require(times == std::vector<double>({0, 1}));
  const auto previous = read();
  require(!imyann::save_trajectory(target.string(), {1, 0}, {1e4, 1e4}, {2e8, 2e8}));
  require(read() == previous);
  fs::remove_all(directory);
}
