python3 -c "
import sqlite3, os
db = '/media/geuba03p/GEURTEN01/cow.db'
print(f'DB size: {os.path.getsize(db)/1e9:.1f} GB')
con = sqlite3.connect(db)
cur = con.cursor()

# Check indexes
cur.execute(\"SELECT name FROM sqlite_master WHERE type='index' AND name LIKE 'idx_%'\")
indexes = [r[0] for r in cur.fetchall()]
print(f'\nIndexes found: {len(indexes)}')
for name in sorted(indexes):
    print(f'  {name}')

# Quick row count on the big table
cur.execute('SELECT COUNT(*) FROM smaxtec_derived')
print(f'\nsmaxtec_derived: {cur.fetchone()[0]:,} rows')
con.close()
"
