#pragma once
#ifndef GENERAL_COMMON_H_
#define GENERAL_COMMON_H_

#ifdef _MSC_VER
#define CONSTEXPR_OR_CONST const
#define CONSTEXPR
#else
#define CONSTEXPR_OR_CONST constexpr
#define CONSTEXPR constexpr
#endif

#endif //GENERAL_COMMON_H