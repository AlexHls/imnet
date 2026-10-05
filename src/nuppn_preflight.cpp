#include "nuppn_preflight.h"
#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <spawn.h>
#include <stdexcept>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <vector>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

extern char **environ;

namespace imyann {
std::filesystem::path nuppn_probe_path() {
#ifdef __APPLE__
  uint32_t size = 0;
  _NSGetExecutablePath(nullptr, &size);
  std::vector<char> buffer(size);
  if (_NSGetExecutablePath(buffer.data(), &size) != 0)
    throw std::runtime_error("Cannot locate the NuPPN startup checker");
  const auto executable = std::filesystem::canonical(buffer.data());
#else
  const auto executable = std::filesystem::canonical("/proc/self/exe");
#endif
  return executable.parent_path() / "imnet_nuppn_probe";
}

void check_nuppn_startup(const std::filesystem::path &run_dir,
                         const std::filesystem::path &probe,
                         std::chrono::milliseconds timeout) {
  namespace fs = std::filesystem;
  std::string pattern = (fs::temp_directory_path() / "imnet-nuppn-check-XXXXXX").string();
  if (!mkdtemp(pattern.data())) throw std::runtime_error("Cannot create NuPPN startup check directory");
  struct Cleanup {
    fs::path path;
    ~Cleanup() { std::error_code ignored; fs::remove_all(path, ignored); }
  } cleanup{pattern};
  const std::string marker = (cleanup.path / "ready").string();
  const std::string log = (cleanup.path / "output.log").string();
  const std::string program = fs::absolute(probe).string();
  const std::string directory = fs::absolute(run_dir).string();
  char *args[] = {const_cast<char *>(program.c_str()), const_cast<char *>(directory.c_str()),
                  const_cast<char *>(marker.c_str()), nullptr};
  posix_spawn_file_actions_t actions;
  int code = posix_spawn_file_actions_init(&actions);
  if (code) throw std::runtime_error("Cannot prepare NuPPN startup check: " + std::string(std::strerror(code)));
  code = posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, log.c_str(),
                                         O_WRONLY | O_CREAT | O_TRUNC, 0600);
  if (!code) code = posix_spawn_file_actions_adddup2(&actions, STDOUT_FILENO, STDERR_FILENO);
  if (!code) code = posix_spawn_file_actions_addopen(&actions, STDIN_FILENO, "/dev/null", O_RDONLY, 0);
  pid_t pid = -1;
  if (!code) code = posix_spawn(&pid, program.c_str(), &actions, nullptr, args, environ);
  posix_spawn_file_actions_destroy(&actions);
  if (code) throw std::runtime_error("Cannot start NuPPN startup checker " + program + ": " + std::strerror(code));

  const auto deadline = std::chrono::steady_clock::now() + timeout;
  int status = 0;
  bool timed_out = false;
  for (;;) {
    const pid_t result = waitpid(pid, &status, WNOHANG);
    if (result == pid) break;
    if (result < 0 && errno != EINTR)
      throw std::runtime_error("Cannot collect NuPPN startup check result");
    if (std::chrono::steady_clock::now() >= deadline) {
      timed_out = true;
      kill(pid, SIGKILL);
      while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
      break;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  std::ifstream ready(marker);
  std::string message;
  std::getline(ready, message);
  if (!timed_out && WIFEXITED(status) && WEXITSTATUS(status) == 0 && message == "imnet-nuppn-ready") return;

  std::ifstream output(log, std::ios::binary | std::ios::ate);
  const auto length = output.tellg();
  if (length > std::streampos(4096)) output.seekg(length - std::streamoff(4096));
  else output.seekg(0);
  const std::string details((std::istreambuf_iterator<char>(output)), {});
  throw std::runtime_error(std::string("NuPPN startup check failed") +
      (timed_out ? " (timed out)" : WIFSIGNALED(status) ? " (native process crashed)" :
       " (native initialization did not complete)") + ":\n" + details);
}
} // namespace imyann
