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
  const auto group = directory / "group";
  fs::create_directory(group);
  const auto first = group / "first", second = group / "second", third = group / "third";
  auto contents = [](const fs::path &path) {
    std::ifstream in(path);
    return std::string(std::istreambuf_iterator<char>(in), {});
  };
  auto writes = [&] {
    return std::vector<imyann::FileWrite>{
        {first.string(), [](auto &out) { out << "new first"; }},
        {second.string(), [](auto &out) { out << "new second"; }},
        {third.string(), [](auto &out) { out << "new third"; }}};
  };
  std::string error;
  require(write_file_atomic(first.string(), [](auto &out) { out << "old first"; }));
  fs::permissions(first, fs::perms::owner_read | fs::perms::owner_write);
  const auto original_permissions = fs::status(first).permissions();
  // Late staging failures must not touch any destination.
  auto batch = writes();
  batch.back().write = [](auto &out) { out << "partial"; out.setstate(std::ios::badbit); };
  require(!imyann::write_files_transactional(batch, error) && !error.empty());
  require(contents(first) == "old first" && !fs::exists(second) && !fs::exists(third));
  require(std::distance(fs::directory_iterator(group), fs::directory_iterator{}) == 1);
  // Force the third rename to fail after two replacements succeeded. Restore
  // the first file and remove the newly-created second file during rollback.
  batch = writes();
  batch.back().write = [&](auto &out) { out << "third"; fs::create_directory(third); };
  require(!imyann::write_files_transactional(batch, error) && !error.empty());
  require(contents(first) == "old first" && !fs::exists(second) && fs::is_directory(third));
  require(fs::status(first).permissions() == original_permissions);
  require(std::distance(fs::directory_iterator(group), fs::directory_iterator{}) == 2);
  fs::remove(third);
  require(imyann::write_files_transactional(writes(), error) && error.empty());
  require(contents(first) == "new first" && contents(second) == "new second" &&
          contents(third) == "new third");
  require(fs::status(first).permissions() == original_permissions);
  require(std::distance(fs::directory_iterator(group), fs::directory_iterator{}) == 3);
  // Aliased destinations cannot accidentally overwrite one another in a group.
  const auto alias = group / "alias";
  fs::create_symlink(first, alias);
  batch = writes();
  batch.back().filename = alias.string();
  require(!imyann::write_files_transactional(batch, error));
  require(contents(first) == "new first" && fs::is_symlink(alias));
  batch = writes();
  batch.front().filename = alias.string();
  require(imyann::write_files_transactional(batch, error));
  require(fs::is_symlink(alias));
  fs::remove_all(directory);
}
