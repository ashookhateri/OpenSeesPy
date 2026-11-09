# راهنمای انتشار openseespyandroid

## ۱. ساخت مخزن در GitHub
1. یک مخزن عمومی (Public) با نام `openseespyandroid` بساز. (runner های arm64 برای مخزن عمومی رایگان هستند.)
2. همه‌ی فایل‌های این پوشه را در آن push کن:
```
git init && git add . && git commit -m "openseespyandroid"
git branch -M main
git remote add origin https://github.com/USERNAME/openseespyandroid.git
git push -u origin main
```

## ۲. تست ساخت (بدون انتشار)
در GitHub: Actions → Build & release openseespyandroid → Run workflow
بعد از اتمام، فایل wheel در بخش Artifacts (dist) هست. می‌توانی آن را روی گوشی دانلود و با
`pip install openseespyandroid-*.whl` تست کنی و بعد `tests/smoke_test.py` را اجرا کنی.

## ۳. اتصال به PyPI (بدون رمز و توکن – Trusted Publishing)
1. در pypi.org حساب بساز / وارد شو.
2. Account settings → Publishing → Add a new pending publisher:
   - PyPI Project Name: `openseespyandroid`
   - Owner: نام کاربری GitHub تو
   - Repository: `openseespyandroid`
   - Workflow: `build-release.yml`
   - Environment: `pypi`
3. در GitHub: Settings → Environments → New environment با نام `pypi`.

## ۴. انتشار نسخه
```
git tag v3.7.1.0
git push --tags
```
Workflow خودش می‌سازد، Release در GitHub می‌سازد و روی PyPI منتشر می‌کند.
برای نسخه‌ی بعد، `version` را در `pyproject.toml` بالا ببر و تگ جدید بزن.

## ۵. نصب روی Pydroid 3
منوی Pip → نوشتن `openseespyandroid` → Install
یا در ترمینال Pydroid: `pip install openseespyandroid`

## نکته‌ها
- اولین ساخت حدود ۱ تا ۳ ساعت طول می‌کشد (کامپایل کل OpenSees برای ۳ نسخه پایتون).
- اگر Pydroid پیغام «No OpenSees build for Python X.Y» داد، نسخه را به `PY_VERSIONS` در workflow اضافه کن.
