/*
    Copyright 2026 Thibaut Paumard

    This file is part of Gyoto.

    Gyoto is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gyoto is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file GyotoAnsiScope.h
 * \brief Color management for terminal output
 *
 * This header provides a complete set of macros for ANSI color and text
 * formatting, designed for consistent, readable terminal output in Gyoto.
 *
 * \section Usage
 *
 * Use these macros with std::cout, std::cerr, or any std::ostream:
 * \code
 * std::cout << GYOTO_ANSI_FG_RED << "Error: "
 *           << GYOTO_ANSI_RESET << "Something went wrong\n";
 * std::cerr << GYOTO_ANSI_BOLD << GYOTO_ANSI_FG_YELLOW << "Warning!"
 *           << GYOTO_ANSI_RESET << "\n";
 * \endcode
 *
 * For automatic color scoping, use with Gyoto::AnsiScope:
 * \code
 * Gyoto::AnsiScope::cout() << GYOTO_ANSI_FG_GREEN << "Success!" << std::endl;
 * \endcode
 *
 * \section Color_Naming
 *
 * The naming scheme is as follows:
 * - GYOTO_ANSI_FG_* : Standard foreground colors (30-37)
 * - GYOTO_ANSI_FG_BRIGHT_* : Bright foreground colors (90-97)
 * - GYOTO_ANSI_BG_* : Standard background colors (40-47)
 * - GYOTO_ANSI_BG_BRIGHT_* : Bright background colors (100-107)
 * - GYOTO_ANSI_* : Text attributes (bold, reverse, blink, etc.)
 * - GYOTO_ANSI_*_OFF : Attribute reset macros
 * - GYOTO_ANSI_RESET : Reset all formatting to default
 *
 * \section Compatibility
 *
 * These macros use standard ANSI escape codes and work on:
 * - Linux terminals (GNOME, Konsole, xterm, etc.)
 * - macOS Terminal and iTerm2
 * - Windows Terminal (since Windows 10 1903)
 * - Most modern terminal emulators
 *
 * For terminals without ANSI support, the codes are ignored (harmless).
 */

#ifndef __GyotoAnsi_H_
#define __GyotoAnsi_H_

#include "GyotoDefs.h"

#include <ostream>

// ============================================================================
// RESET TO DEFAULT
// ============================================================================

/// \brief Reset all attributes (color, bold, underline, etc.) to default
#define GYOTO_ANSI_RESET "\033[0m"

// ============================================================================
// FOREGROUND COLORS (Standard: 30-37)
// ============================================================================

/// \brief Standard foreground: Black
#define GYOTO_ANSI_FG_BLACK "\033[30m"
/// \brief Standard foreground: Red
#define GYOTO_ANSI_FG_RED "\033[31m"
/// \brief Standard foreground: Green
#define GYOTO_ANSI_FG_GREEN "\033[32m"
/// \brief Standard foreground: Yellow
#define GYOTO_ANSI_FG_YELLOW "\033[33m"
/// \brief Standard foreground: Blue
#define GYOTO_ANSI_FG_BLUE "\033[34m"
/// \brief Standard foreground: Magenta
#define GYOTO_ANSI_FG_MAGENTA "\033[35m"
/// \brief Standard foreground: Cyan
#define GYOTO_ANSI_FG_CYAN "\033[36m"
/// \brief Standard foreground: White
#define GYOTO_ANSI_FG_WHITE "\033[37m"

// ============================================================================
// FOREGROUND COLORS (Bright: 90-97)
// ============================================================================

/// \brief Bright foreground: Black (dark gray on some terminals)
#define GYOTO_ANSI_FG_BRIGHT_BLACK "\033[90m"
/// \brief Bright foreground: Red
#define GYOTO_ANSI_FG_BRIGHT_RED "\033[91m"
/// \brief Bright foreground: Green
#define GYOTO_ANSI_FG_BRIGHT_GREEN "\033[92m"
/// \brief Bright foreground: Yellow
#define GYOTO_ANSI_FG_BRIGHT_YELLOW "\033[93m"
/// \brief Bright foreground: Blue
#define GYOTO_ANSI_FG_BRIGHT_BLUE "\033[94m"
/// \brief Bright foreground: Magenta
#define GYOTO_ANSI_FG_BRIGHT_MAGENTA "\033[95m"
/// \brief Bright foreground: Cyan
#define GYOTO_ANSI_FG_BRIGHT_CYAN "\033[96m"
/// \brief Bright foreground: White
#define GYOTO_ANSI_FG_BRIGHT_WHITE "\033[97m"

// ============================================================================
// BACKGROUND COLORS (Standard: 40-47)
// ============================================================================

/// \brief Standard background: Black
#define GYOTO_ANSI_BG_BLACK "\033[40m"
/// \brief Standard background: Red
#define GYOTO_ANSI_BG_RED "\033[41m"
/// \brief Standard background: Green
#define GYOTO_ANSI_BG_GREEN "\033[42m"
/// \brief Standard background: Yellow
#define GYOTO_ANSI_BG_YELLOW "\033[43m"
/// \brief Standard background: Blue
#define GYOTO_ANSI_BG_BLUE "\033[44m"
/// \brief Standard background: Magenta
#define GYOTO_ANSI_BG_MAGENTA "\033[45m"
/// \brief Standard background: Cyan
#define GYOTO_ANSI_BG_CYAN "\033[46m"
/// \brief Standard background: White
#define GYOTO_ANSI_BG_WHITE "\033[47m"

// ============================================================================
// BACKGROUND COLORS (Bright: 100-107)
// ============================================================================

/// \brief Bright background: Black (dark gray on some terminals)
#define GYOTO_ANSI_BG_BRIGHT_BLACK "\033[100m"
/// \brief Bright background: Red
#define GYOTO_ANSI_BG_BRIGHT_RED "\033[101m"
/// \brief Bright background: Green
#define GYOTO_ANSI_BG_BRIGHT_GREEN "\033[102m"
/// \brief Bright background: Yellow
#define GYOTO_ANSI_BG_BRIGHT_YELLOW "\033[103m"
/// \brief Bright background: Blue
#define GYOTO_ANSI_BG_BRIGHT_BLUE "\033[104m"
/// \brief Bright background: Magenta
#define GYOTO_ANSI_BG_BRIGHT_MAGENTA "\033[105m"
/// \brief Bright background: Cyan
#define GYOTO_ANSI_BG_BRIGHT_CYAN "\033[106m"
/// \brief Bright background: White
#define GYOTO_ANSI_BG_BRIGHT_WHITE "\033[107m"

// ============================================================================
// TEXT ATTRIBUTES
// ============================================================================

/// \brief Bold text
#define GYOTO_ANSI_BOLD "\033[1m"
/// \brief Disable bold
#define GYOTO_ANSI_BOLD_OFF "\033[22m"

/// \brief Dim/faint text
#define GYOTO_ANSI_DIM "\033[2m"
/// \brief Disable dim/faint
#define GYOTO_ANSI_DIM_OFF "\033[22m"

/// \brief Italic text
#define GYOTO_ANSI_ITALIC "\033[3m"
/// \brief Disable italic
#define GYOTO_ANSI_ITALIC_OFF "\033[23m"

/// \brief Underlined text
#define GYOTO_ANSI_UNDERLINE "\033[4m"
/// \brief Disable underline
#define GYOTO_ANSI_UNDERLINE_OFF "\033[24m"

/// \brief Blinking text
#define GYOTO_ANSI_BLINK "\033[5m"
/// \brief Disable blink
#define GYOTO_ANSI_BLINK_OFF "\033[25m"

/// \brief Reverse video (swap foreground and background)
#define GYOTO_ANSI_REVERSE "\033[7m"
/// \brief Disable reverse video
#define GYOTO_ANSI_REVERSE_OFF "\033[27m"

/// \brief Hidden text (useful for passwords)
#define GYOTO_ANSI_HIDDEN "\033[8m"
/// \brief Show hidden text
#define GYOTO_ANSI_HIDDEN_OFF "\033[28m"

/// \brief Strikethrough text
#define GYOTO_ANSI_STRIKETHROUGH "\033[9m"
/// \brief Disable strikethrough
#define GYOTO_ANSI_STRIKETHROUGH_OFF "\033[29m"

// ============================================================================
// COMBINED FORMATTING (Common use cases)
// ============================================================================

/// \brief Format "SEVERE:" tag
#define GYOTO_ANSI_SEVERE_TAG					\
  GYOTO_ANSI_RESET GYOTO_ANSI_BOLD GYOTO_ANSI_FG_BRIGHT_RED
/// \brief Format severe warning
#define GYOTO_ANSI_SEVERE			\
  GYOTO_ANSI_RESET GYOTO_ANSI_FG_RED
/// \brief Format "WARNING:" tag
#define GYOTO_ANSI_WARNING_TAG \
  GYOTO_ANSI_RESET GYOTO_ANSI_BOLD GYOTO_ANSI_FG_BRIGHT_YELLOW
/// \brief Format warning message
#define GYOTO_ANSI_WARNING GYOTO_ANSI_RESET GYOTO_ANSI_FG_YELLOW
/// \brief Format "INFO:" tag
#define GYOTO_ANSI_INFO_TAG					\
  GYOTO_ANSI_RESET GYOTO_ANSI_BOLD GYOTO_ANSI_FG_BRIGHT_GREEN
/// \brief Format info message
#define GYOTO_ANSI_INFO GYOTO_ANSI_RESET GYOTO_ANSI_FG_GREEN

/// \brief Format "DEBUG:" tag
#define GYOTO_ANSI_DEBUG_TAG					\
  GYOTO_ANSI_RESET GYOTO_ANSI_BOLD GYOTO_ANSI_FG_BRIGHT_BLUE

/// \brief Format pretty function in debug message
#define GYOTO_ANSI_DEBUG_PRETTY_FUNCTION	\
  GYOTO_ANSI_RESET GYOTO_ANSI_FG_BRIGHT_CYAN

/// \brief Format debug message
#define GYOTO_ANSI_DEBUG			\
  GYOTO_ANSI_RESET GYOTO_ANSI_FG_BLUE

/// \brief Format value in debug message
#define GYOTO_ANSI_DEBUG_VALUE			\
  GYOTO_ANSI_RESET GYOTO_ANSI_FG_BRIGHT_YELLOW

/// \brief Format variable in debug message
#define GYOTO_ANSI_DEBUG_VARIABLE			\
  GYOTO_ANSI_RESET GYOTO_ANSI_FG_BRIGHT_MAGENTA

namespace Gyoto {
  /**
   * \brief Scoped color management for terminal output
   *
   * This class sets a color on construction and resets to default
   * on destruction, ensuring the stream is always left in a clean state.
   */
  class AnsiScope;

}

class Gyoto::AnsiScope {
  std::ostream& stream_;
  bool reset_;
public:
  AnsiScope(std::ostream& s);
  ~AnsiScope();
  void no_reset();
  static Gyoto::AnsiScope cout();
  static Gyoto::AnsiScope cerr();
  void append(std::string txt);
  template<typename T>
  AnsiScope& operator<<(const T& value) {
    stream_ << value;
    return *this;
  }
  AnsiScope& operator<<(std::ostream& (*manip)(std::ostream&));
  void flush();
};
#endif
