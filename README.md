## 專案簡介
本專案為 **MMIP FinalProject**：實作 **IR / VIS 影像融合**，提供三種分解方法並以相同的融合規則產生融合結果：

- **DWT**（Discrete Wavelet Transform）
- **SWT**（Stationary Wavelet Transform / Undecimated）
- **DT-CWT**（Dual-Tree Complex Wavelet Transform, 2D）

融合流程（各模式一致）：
- **VIS RGB → YCbCr**，只對 **Y** 做分解/融合
- **IR 灰階 → double**，做分解/融合
- 融合後得到 **Y_fused**，再與原本 **Cb/Cr** 合成回 RGB
- 輸出 `PNG` 到使用者指定資料夾

> 融合子頻規則：低頻/高頻皆採用 **local energy weighting**（3×3 window）。

---

## 專案結構
```
FinalProject/
  include/
    dwt.h        # DWT 入口宣告
    swt.h        # SWT 入口宣告
    dtcwt.h      # DT-CWT 資料結構/入口宣告
    utils.h      # 共用 I/O、色彩轉換、批次掃檔等
  src/
    main.cpp     # CLI 參數解析 + 模式分派
    dwt.cpp      # DWT 融合
    swt.cpp      # SWT 融合
    dtcwt.cpp    # DT-CWT forward/inverse + 融合 + PR test
    utils.cpp    # stb_image/stb_image_write 實作 + 共用工具
  external/
    wavelib/     # 2D DWT/SWT 使用的第三方庫
  assets/
    Vis/         # 可選：VIS 資料集（png）
    Ir/          # 可選：IR  資料集（png）
```

---

## 需求環境
- **Windows 10/11**
- **Visual Studio 2022**（MSVC）
- **CMake ≥ 3.20**
- C++17

## 第三方套件 / 相依
- **wavelib**：本專案的 `DWT` / `SWT` 由 `wavelib` 提供（來源：[rafat/wavelib](https://github.com/rafat/wavelib)），專案內以 `external/wavelib` 方式引入並由 CMake `add_subdirectory` 編譯/連結。
- **stb_image / stb_image_write**：用於 PNG 影像的讀取與寫入。
  標頭檔放置於 `include/`，並於 `src/utils.cpp` 中定義對應的 `STB_*_IMPLEMENTATION` 以產生實體實作。

---

## 編譯方式（Windows / Visual Studio）
在專案根目錄執行：

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

執行檔位置：
- `build\Release\FinalProject.exe`

---

## 執行方式（CLI）
所有模式都要求：
- 單張：提供 `VIS` 與 `IR` 的檔案路徑
- 批次（DWT/SWT/DT-CWT）：提供 `VIS` 資料夾與 `IR` 資料夾
- **必須指定** `--output <output_dir>`

### 單張融合
```powershell
build\Release\FinalProject.exe --mode dwt   assets/Vis/00955.png assets/Ir/00955.png --output build/out_dwt
build\Release\FinalProject.exe --mode swt   assets/Vis/00955.png assets/Ir/00955.png --output build/out_swt
build\Release\FinalProject.exe --mode dtcwt assets/Vis/00955.png assets/Ir/00955.png --output build/out_dtcwt
```

輸出檔名格式（各 mode 一致）：
- `fused_<vis_filename>_<mode>.png`

### 批次融合
```powershell
build\Release\FinalProject.exe --mode dwt --batch assets/Vis assets/Ir --output build/out_dwt
build\Release\FinalProject.exe --mode swt --batch assets/Vis assets/Ir --output build/out_swt
build\Release\FinalProject.exe --mode dtcwt --batch assets/Vis assets/Ir --output build/out_dtcwt
```

---

## 參數說明
- **`--mode <dwt|swt|dtcwt|dtcwt_test>`**
- **`--batch`**：啟用批次模式
- **`--output <dir>`**：輸出資料夾（必填）
---



