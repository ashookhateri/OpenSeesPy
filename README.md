# OpenSeesPy for Android — `openseespyandroid`

[English](#english) | [فارسی](#فارسی)

[![PyPI](https://img.shields.io/pypi/v/openseespyandroid)](https://pypi.org/project/openseespyandroid/)
![Platform](https://img.shields.io/badge/platform-Android%20arm64-green)
![Python](https://img.shields.io/badge/python-3.10%20to%203.14-blue)

---

## English

### What is this?

**openseespyandroid** brings [OpenSees](https://github.com/OpenSees/OpenSees) — the Open System for Earthquake Engineering Simulation — to Android phones and tablets as a ready-to-install Python package.

OpenSees is a finite element framework used by structural and geotechnical engineers for linear and nonlinear static, dynamic, pushover and time-history analysis. The official [OpenSeesPy](https://github.com/zhuminjie/OpenSeesPy) package only ships builds for Windows, macOS and Linux x86_64, so it cannot run on a phone. This project compiles OpenSees natively for **Android arm64 (aarch64)** and packages it so that a single `pip install` works inside Android Python apps — no root, no PC, no Linux emulator.

### Author

**Dr. Amir Khateri** — GitHub: [@ashookhateri](https://github.com/ashookhateri)

### Features

- Same API as `openseespy.opensees`: existing scripts work after changing one import line
- Native arm64 build (C, C++ and Fortran compiled for Android bionic), full speed on the device
- One wheel for Python **3.10, 3.11, 3.12, 3.13 and 3.14**; the loader picks the right build automatically
- All required shared libraries bundled inside the package; nothing else to install
- Runs fully offline once installed

### Android apps you can use it with

| App | Where to get it | Notes |
|---|---|---|
| **Pydroid 3** (recommended) | [Google Play](https://play.google.com/store/apps/details?id=ru.iiec.pydroid3) | Easiest. Install from the app's Pip menu or its terminal. Includes an editor, console and Matplotlib for plotting results. |
| **Termux** | [F-Droid](https://f-droid.org/packages/com.termux/) or [GitHub](https://github.com/termux/termux-app) | Full Linux-style terminal. Run `pkg install python` first. Best for long analyses and batch scripts. Use the F-Droid/GitHub build, not the old Google Play one. |
| **VSCodroid** | [GitHub](https://github.com/rmyndharis/VSCodroid) / [Google Play](https://play.google.com/store/apps/details?id=com.vscodroid) | Full VS Code on Android with a real terminal and Python 3.14 + pip built in. Run `pip install openseespyandroid` in its terminal and use the Run button. Needs Android 11+. |
| **QPython 3L** | Google Play | Works when its bundled Python is 3.10 or newer and pip is available. |
| **Code editors with Termux** (e.g. Acode) | Google Play / F-Droid | Write scripts in the editor and run them in Termux. |

Requirements: Android 7.0 (API 24) or newer, a 64-bit ARM (arm64-v8a) device. Almost every phone from the last several years qualifies. 32-bit ARM and x86 devices are not supported.

### Install

**Pydroid 3:** open the menu → **Pip** → type `openseespyandroid` → **Install**.
If the install fails with a cache error, open **Terminal** in Pydroid and run:

```
pip install --no-cache-dir openseespyandroid
```

**Termux:**

```
pkg update && pkg install python
pip install openseespyandroid
```

**Upgrade to the latest version:**

```
pip install --no-cache-dir --upgrade openseespyandroid
```

### Use

```python
import openseespyandroid.opensees as ops

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)
ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 3.0)
ops.fix(1, 1, 1, 1)
ops.geomTransf('Linear', 1)
ops.element('elasticBeamColumn', 1, 1, 2, 0.09, 2.0e8, 6.75e-4, 1)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 10.0, 0.0, 0.0)
ops.system('BandGeneral'); ops.numberer('RCM'); ops.constraints('Plain')
ops.integrator('LoadControl', 1.0); ops.algorithm('Linear'); ops.analysis('Static')
ops.analyze(1)
print('tip displacement =', ops.nodeDisp(2, 1))
```

**Porting an existing OpenSeesPy script:** replace

```python
import openseespy.opensees as ops
```

with

```python
import openseespyandroid.opensees as ops
```

Everything else stays the same. A self-checking test is in [`tests/smoke_test.py`](tests/smoke_test.py).

### Troubleshooting

| Message | Fix |
|---|---|
| `No module named 'openseespy.opensees'` | You are using the desktop import. Use `import openseespyandroid.opensees as ops`. |
| `RuntimeError` saying the platform must be x86_64 | The desktop `openseespy` package is installed. Run `pip uninstall openseespy openseespylinux`, then install `openseespyandroid`. |
| `library "libpython3.X.so" not found` | An old version is installed. Run `pip uninstall -y openseespyandroid`, then `pip install --no-cache-dir openseespyandroid`, and restart the app. |
| `No OpenSees build for Python 3.X` | Your app's Python version is not built yet. Open an issue with the version number. |
| `OSError: [Errno 2]` during install in Pydroid | Pydroid cache problem. Install from Pydroid's terminal with `--no-cache-dir`. |

### How it is built

The GitHub Actions workflow [`.github/workflows/build-release.yml`](.github/workflows/build-release.yml) runs on a native arm64 runner inside the official Termux Docker image (a real Android bionic userland). It compiles OpenSees with clang and LLVM flang for each Python version, links each module only against its own `libpython`, bundles dependencies with an `$ORIGIN` RUNPATH, and publishes a `py3-none-any` wheel to PyPI and GitHub Releases whenever a `v*` tag is pushed. A Persian step-by-step build guide is in [`BUILD_GUIDE_FA.md`](BUILD_GUIDE_FA.md).

### License

OpenSees is © The Regents of the University of California. See [`LICENSE`](LICENSE). This package redistributes compiled binaries of the upstream source with Android build scripts. It is a community project and not affiliated with the official OpenSees or OpenSeesPy teams.

---

<div dir="rtl">

## فارسی

### این پروژه چیست؟

**openseespyandroid** نرم‌افزار [OpenSees](https://github.com/OpenSees/OpenSees) (سامانه‌ی باز شبیه‌سازی مهندسی زلزله) را به‌صورت یک بسته‌ی پایتون آماده‌ی نصب، روی گوشی و تبلت اندرویدی می‌آورد.

OpenSees یک چارچوب اجزای محدود است که مهندسان سازه و ژئوتکنیک برای تحلیل‌های خطی و غیرخطی، استاتیکی، دینامیکی، پوش‌آور و تاریخچه‌ی زمانی از آن استفاده می‌کنند. بسته‌ی رسمی [OpenSeesPy](https://github.com/zhuminjie/OpenSeesPy) فقط برای ویندوز، مک و لینوکس x86_64 ساخته می‌شود و روی گوشی اجرا نمی‌شود. این پروژه OpenSees را به‌صورت بومی برای **اندروید arm64 (aarch64)** کامپایل می‌کند تا فقط با یک دستور `pip install` داخل برنامه‌های پایتون اندروید کار کند؛ بدون روت، بدون کامپیوتر و بدون شبیه‌ساز لینوکس.

### سازنده

**دکتر امیر خاطری (Dr. Amir Khateri)** — گیت‌هاب: [@ashookhateri](https://github.com/ashookhateri)

### ویژگی‌ها

- همان API بسته‌ی `openseespy.opensees`؛ اسکریپت‌های قبلی فقط با تغییر یک خط import اجرا می‌شوند
- ساخت بومی arm64 (کدهای C، ‏C++ و فرترن برای اندروید)، با سرعت کامل روی خود گوشی
- یک بسته برای پایتون **3.10، 3.11، 3.12، 3.13 و 3.14**؛ نسخه‌ی مناسب خودکار انتخاب می‌شود
- همه‌ی کتابخانه‌های لازم داخل بسته هستند و چیز دیگری لازم نیست نصب شود
- بعد از نصب، کاملاً آفلاین کار می‌کند

### برنامه‌های اندرویدی قابل استفاده

| برنامه | محل دریافت | توضیح |
|---|---|---|
| **Pydroid 3** (پیشنهادی) | [گوگل‌پلی](https://play.google.com/store/apps/details?id=ru.iiec.pydroid3) | ساده‌ترین راه. نصب از منوی Pip یا ترمینال برنامه. ویرایشگر، کنسول و Matplotlib برای رسم نمودار نتایج دارد. |
| **Termux** | [F-Droid](https://f-droid.org/packages/com.termux/) یا [گیت‌هاب](https://github.com/termux/termux-app) | ترمینال کامل شبیه لینوکس. اول `pkg install python` را اجرا کنید. برای تحلیل‌های طولانی و اجرای دسته‌ای بهترین گزینه است. نسخه‌ی F-Droid یا گیت‌هاب را نصب کنید، نه نسخه‌ی قدیمی گوگل‌پلی. |
| **VSCodroid** | [گیت‌هاب](https://github.com/rmyndharis/VSCodroid) / [گوگل‌پلی](https://play.google.com/store/apps/details?id=com.vscodroid) | خود VS Code روی اندروید با ترمینال واقعی و پایتون 3.14 و pip داخلی. در ترمینالش `pip install openseespyandroid` را بزنید و با دکمه‌ی Run اجرا کنید. اندروید 11 به بالا. |
| **QPython 3L** | گوگل‌پلی | اگر پایتون داخلش 3.10 یا جدیدتر باشد و pip داشته باشد کار می‌کند. |
| **ویرایشگرهای کد همراه Termux** (مثل Acode) | گوگل‌پلی / F-Droid | اسکریپت را در ویرایشگر بنویسید و در Termux اجرا کنید. |

پیش‌نیاز: اندروید 7.0 (API 24) یا بالاتر و گوشی 64 بیتی ARM ‏(arm64-v8a). تقریباً همه‌ی گوشی‌های چند سال اخیر این شرایط را دارند. دستگاه‌های ARM ‏32 بیتی و x86 پشتیبانی نمی‌شوند.

### نصب

**در Pydroid 3:** منو ← **Pip** ← نوشتن `openseespyandroid` ← **Install**.
اگر نصب با خطای cache متوقف شد، **Terminal** را در Pydroid باز کنید و بزنید:

</div>

```
pip install --no-cache-dir openseespyandroid
```

<div dir="rtl">

**در Termux:**

</div>

```
pkg update && pkg install python
pip install openseespyandroid
```

<div dir="rtl">

**به‌روزرسانی به آخرین نسخه:**

</div>

```
pip install --no-cache-dir --upgrade openseespyandroid
```

<div dir="rtl">

### استفاده

</div>

```python
import openseespyandroid.opensees as ops

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)
# ... بقیه‌ی دستورها دقیقاً مثل openseespy
```

<div dir="rtl">

**تبدیل اسکریپت‌های قبلی OpenSeesPy:** فقط خط

</div>

```python
import openseespy.opensees as ops
```

<div dir="rtl">

را با این خط عوض کنید:

</div>

```python
import openseespyandroid.opensees as ops
```

<div dir="rtl">

بقیه‌ی کد هیچ تغییری نمی‌خواهد. یک تست آماده در [`tests/smoke_test.py`](tests/smoke_test.py) هست.

### رفع اشکال

| پیام خطا | راه‌حل |
|---|---|
| `No module named 'openseespy.opensees'` | از import نسخه‌ی دسکتاپ استفاده کرده‌اید. بنویسید `import openseespyandroid.opensees as ops` |
| خطای `RuntimeError` درباره‌ی x86_64 | بسته‌ی دسکتاپ `openseespy` نصب است. اول `pip uninstall openseespy openseespylinux` و بعد `openseespyandroid` را نصب کنید. |
| `library "libpython3.X.so" not found` | نسخه‌ی قدیمی نصب است. `pip uninstall -y openseespyandroid` و بعد `pip install --no-cache-dir openseespyandroid` را بزنید و برنامه را دوباره باز کنید. |
| `No OpenSees build for Python 3.X` | نسخه‌ی پایتون برنامه‌ی شما هنوز ساخته نشده. با ذکر شماره‌ی نسخه یک Issue باز کنید. |
| `OSError: [Errno 2]` هنگام نصب در Pydroid | مشکل cache خود Pydroid است. از ترمینال Pydroid با `--no-cache-dir` نصب کنید. |

### روش ساخت

فایل [`.github/workflows/build-release.yml`](.github/workflows/build-release.yml) روی رانر بومی arm64 گیت‌هاب و داخل ایمیج رسمی Termux (محیط واقعی اندروید) اجرا می‌شود؛ OpenSees را با clang و LLVM flang برای هر نسخه‌ی پایتون می‌سازد، هر ماژول را فقط به `libpython` نسخه‌ی خودش لینک می‌کند و با هر تگ `v*` بسته را در PyPI و GitHub Releases منتشر می‌کند. راهنمای کامل ساخت به فارسی در [`BUILD_GUIDE_FA.md`](BUILD_GUIDE_FA.md) است.

### مجوز

حق نشر OpenSees متعلق به The Regents of the University of California است؛ فایل [`LICENSE`](LICENSE) را ببینید. این بسته فایل‌های کامپایل‌شده‌ی کد اصلی را همراه اسکریپت‌های ساخت اندروید منتشر می‌کند. این یک پروژه‌ی مستقل است و وابسته به تیم رسمی OpenSees یا OpenSeesPy نیست.

</div>
