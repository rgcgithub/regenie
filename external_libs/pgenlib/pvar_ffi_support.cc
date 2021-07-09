// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.
//
//
// *  Modified by Joelle Mbatchou - Apr 27 2021
// *  - kept only needed functions and modify to recude needed headers

#include "pvar_ffi_support.h"

namespace plink2 {

RefcountedWptr* CreateRefcountedWptr(uintptr_t size) {
  RefcountedWptr* rwp = static_cast<RefcountedWptr*>(malloc(sizeof(RefcountedWptr)));
  if (!rwp) {
    return nullptr;
  }
  rwp->ref_ct = 1;
  rwp->p = static_cast<uintptr_t*>(malloc(sizeof(uintptr_t) * size));

  return rwp;
}

void CondReleaseRefcountedWptr(RefcountedWptr** rwpp) {
  RefcountedWptr* rwp = *rwpp;
  if (!rwp) {
    return;
  }
  --rwp->ref_ct;
  if (!rwp->ref_ct) {
    free(rwp->p);
    free(rwp);
  }
  *rwpp = nullptr;
}

}  // namespace plink2
