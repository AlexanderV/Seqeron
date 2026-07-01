# Claude Code Skills Strategy — Seqeron

> **Мета.** Перетворити добре документовану біобібліотеку `Seqeron.Genomics` (+`SuffixTree`)
> та 11 MCP-серверів (427 інструментів) на систему, яка **надійно розв'язує складні комплексні
> біологічні задачі** — як через MCP-оркестрацію (LLM-клієнт), так і через прямий C#/.NET API.
>
> **Статус:** DRAFT / plan of record. Реалізація — фазами (див. §8).
> **Останнє оновлення:** 2026-07-01.

---

## 1. Рішення, зафіксовані перед стартом

| Питання | Рішення |
|---|---|
| Аудиторія | **Обидва профілі**: (a) кінцевий користувач через MCP; (b) розробник через C# API. Кожен доменний скіл дає рецепт для обох режимів. |
| Архітектура | **Багато незалежних скілів** (не єдиний hub). Кожен вмикається окремим trigger'ом. |
| Проти дрейфу | Обов'язковий: скіли **посилаються** на source of truth і мають **спільний generator** таблиць (§5). |
| Перший артефакт | Цей `STRATEGY.md`. |

---

## 2. Що вже існує (і чому це визначає підхід)

- **Бібліотека:** 258 алгоритмічних юнітів; ~15k тестів; 10 методологій тестування;
  campaign валідації + `docs/Validation/LIMITATIONS.md` + рантаймовий `LimitationPolicy`
  (Strict < Moderate < Permissive, дефолт Moderate; 9 guarded-юнітів).
- **MCP:** 11 серверів, **427 інструментів**, кожен — gold-standard binding + схема + тест.
  Source of truth статусу: `docs/mcp/MCP_STATUS.md`.
- **Документація (вже вичерпна):**
  - `docs/algorithms/<Area>/<Unit>.md` — 247 доків з інваріантами, формулами, контрактами, `Test Unit ID`.
  - `docs/mcp/tools/<server>/<tool>.md` (+`.mcp.json`) — 427 per-tool доків: схема I/O, `Method ID`, лінк на джерело.
  - `docs/checklists/*.md` — трекери тест-дисциплін.

**Висновок.** Вузьке місце для розв'язання задач — **не документація і не покриття**, а три речі,
які довідник сам по собі не закриває:

1. **Discovery / routing** — з 427 інструментів обрати правильні; підключення всіх 11 серверів
   одночасно роздуває контекст.
2. **Orchestration** — коректний порядок кроків, вибір параметрів, передача форматів
   (FASTA → aligner → stats), 0-based координати, одиниці.
3. **Наукова валідність** — рахувати інструментами, а не «в умі»; не виходити за validated envelope
   (`LIMITATIONS.md` / `LimitationPolicy`); крос-перевіряти; вести provenance.

Скіли — це **тонкий шар маршрутизації + дисципліни поверх source of truth**, а не переписування доків.

---

## 3. Принципи дизайну (обов'язкові для кожного скіла)

1. **Point, don't duplicate.** Скіл дає _рецепт_ і _лінки_ на `docs/mcp/tools/...`, `docs/algorithms/...`,
   `LIMITATIONS.md`. Схеми/формули не копіюємо — вони дрейфують. Максимум — стислий індекс імен інструментів,
   згенерований автоматично (§5).
2. **Self-updating.** Будь-які таблиці «сервер → домен → інструменти» породжуються скриптом з
   `docs/mcp/tools/` + `MCP_STATUS.md`. Скіл із застарілою таблицею = баг збірки.
3. **Dual-mode.** Кожен доменний скіл містить дві доріжки: **MCP** (імена інструментів + порядок викликів)
   та **C# API** (`Method ID` → реальний виклик `Seqeron.Genomics`, з `LimitationPolicy` bootstrap де треба).
4. **Rigor-by-default.** Доменні скіли делегують наскрізні правила скілу `bio-rigor` (§4.1), не повторюють їх.
5. **Progressive disclosure.** `SKILL.md` короткий (trigger + дерево рішень + рецепти-кістяки);
   деталі — у reference-файлах скіла, що вантажаться on-demand.
6. **Envelope-aware.** Якщо задача чіпає guarded-юніт поза його `MinimumMode` — скіл ЗУПИНЯЄ і повідомляє,
   а не «продавлює» результат.
7. **Reproducible output.** Результат супроводжується provenance: перелік інструментів/`Method ID` + параметри.

---

## 4. Каталог скілів (незалежні)

Позначки режиму: **[MCP]** / **[API]** / **[both]**.

### 4.1 Наскрізні (cross-cutting)

- **`bio-rigor`** [both] — guardrail, вмикається за замовчуванням при роботі з даними.
  Зміст: tool-only дисципліна (парсити/рахувати інструментами, не вручну) · provenance-ланцюг ·
  повага до `LimitationPolicy`/`LIMITATIONS.md` · alpha/clinical-disclaimer ·
  крос-звірка критичних результатів двома незалежними шляхами · перевірка одиниць і 0-based координат.
- **`seqeron-discovery`** [both] — «чи є інструмент/алгоритм для X?».
  Скрипт-пошук по `docs/mcp/tools/` + `docs/algorithms/` → повертає ім'я інструмента, `Method ID`,
  лінк на схему й на алгоритм-док. Дешевше, ніж тримати 427 схем у контексті.

### 4.2 Доменні (workflow-родини, а не «сервер=скіл»)

| Скіл | Покриває сервери | Типові задачі |
|---|---|---|
| **`bio-qc`** | Sequence, Parsers | валідація DNA/RNA/білка, GC%, Tm, композиція, складність, парсинг форматів |
| **`bio-alignment`** | Alignment, Core | pairwise (NW/SW/semi-global), MSA, alignment-stats, similarity/edit-distance |
| **`bio-assembly`** | Analysis, Core | k-mer, de Bruijn/OLC, coverage, assembly-stats, повтори |
| **`bio-annotation`** | Annotation, Analysis | гени/ORF/промотори, мотиви, variant calling, ефект варіантів |
| **`bio-phylo-popgen`** | Phylogenetics, Population | дистанції, побудова дерев, phylo-статистика, popgen-метрики, тести відбору |
| **`bio-metagenomics`** | Metagenomics | таксономічна класифікація, профілювання спільнот |
| **`bio-moldesign`** | MolTools | праймери/зонди, CRISPR-гайди, codon-оптимізація, рестрикція |
| **`bio-chromosome`** | Chromosome | хромосомний рівень, структурні варіанти, GC-skew |

Кожен доменний скіл: `SKILL.md` (trigger + коли-вмикати + 3–6 канонічних пайплайнів як кістяки) +
reference-файли (детальні рецепти, gotchas, вибір параметрів) + автозгенерований індекс інструментів свого домену.

### 4.3 Для профілю «розробник»

- **`seqeron-dev`** [API] — конвенції роботи з бібліотекою напряму: `LimitationPolicy` bootstrap для тестів,
  `TryCreate` vs ctor, guarded-юніти та їх `MinimumMode`, посилання на `docs/algorithms` й тест-дисципліни.
  (Дев-скіли для _розвитку_ бібліотеки — тести/архітектура — вже покриті наявними `clean-architecture`,
  `clean-code`, `simian` і кампаніями в MEMORY; сюди не дублюємо.)

---

## 5. Механізм проти дрейфу (критично для «багато незалежних скілів»)

Ризик обраної архітектури — N скілів, що розходяться з реальним набором інструментів. Контрзаходи:

1. **Один generator** `scripts/skills/gen-catalog.*`:
   читає `docs/mcp/tools/**` + `docs/mcp/MCP_STATUS.md` → генерує на кожен домен файл
   `_generated/tools.md` (таблиця: tool · server · Method ID · лінк). Скіли `@include`/посилаються на нього.
2. **Границі генерації** — згенеровані блоки обгорнуті маркерами `<!-- BEGIN generated -->…<!-- END generated -->`;
   ручний текст поза ними не чіпається.
3. **CI-перевірка** — окремий тест/скрипт `check-catalog-fresh`: падає, якщо `_generated/` розійшлось із source.
   Додати в наявний тест-pipeline поряд з MCP-контракт-тестами.
4. **Правило рев'ю** — новий MCP-інструмент вважається «done» лише коли перегенеровано каталог скілів.

---

## 6. Конвенції файлів скілів

```
.claude/skills/<skill-name>/
  SKILL.md            # frontmatter: name, description (trigger!). Короткий: дерево рішень + кістяки.
  reference/*.md      # детальні рецепти, вантажаться on-demand
  _generated/         # автозгенеровані таблиці інструментів (не редагувати вручну)
  scripts/*           # опційні хелпери (discovery-пошук тощо)
```

- `description` у frontmatter — це **trigger**: конкретні дієслова/іменники біозадачі
  («вирівняй», «поклич варіанти», «спроєктуй праймери», «побудуй дерево»), щоб скіл вмикався автоматично.
- SKILL.md ≤ ~150 рядків; усе інше — progressive disclosure.
- Кожен рецепт закінчується блоком **Provenance** (які інструменти/`Method ID`, які параметри).

---

## 7. Критерії «ідеального» скіла (acceptance)

Скіл вважається готовим, коли:

1. **Triggering** — вмикається на реалістичних формулюваннях задачі домену (перевірено ≥3 промптами).
2. **Routing** — обирає правильні інструменти без завантаження всіх 427 схем.
3. **Both-mode** — дає робочий MCP-рецепт і робочий C#-рецепт для тієї ж задачі.
4. **Rigor** — делегує `bio-rigor`; поважає envelope; віддає provenance.
5. **No-drift** — індекс інструментів згенерований, `check-catalog-fresh` зелений.
6. **Grounded** — відкатаний на ≥1 реальному end-to-end workflow (взяти з `docs/mcp/README.md` прикладів
   і розширити до «складної» багатокрокової задачі).

---

## 8. Фази реалізації

- **Ф0 — фундамент. ✅ ГОТОВО (2026-07-01).** `scripts/skills/gen-catalog.py` (+`--check`,
  `docs/skills/domain-map.json`, `docs/skills/_generated/{catalog.json,tool-catalog.md}`),
  `scripts/skills/check-catalog-fresh.sh`; `scripts/skills/find-tool.py`;
  скіли `.claude/skills/bio-rigor/` та `.claude/skills/seqeron-discovery/`.
  427 інструментів; drift-детекція перевірена; обидва скіли зареєстровані harness'ом.
- **Ф1 — найцінніші домени. ✅ ГОТОВО (2026-07-01).** `bio-annotation` (annotation+analysis, 188 tools),
  `bio-moldesign` (moltools, 47), `bio-alignment` (alignment+core, 34) — dual-mode, end-to-end кейси,
  guarded-юніти з STOP-правилом; `_generated/tools.md` slice'и згенеровано (34/47/188), freshness зелений.
- **Ф2 — решта доменів.** `bio-qc`, `bio-assembly`, `bio-phylo-popgen`, `bio-metagenomics`, `bio-chromosome`.
- **Ф3 — розробник + поліш.** `seqeron-dev`; validation-звірки; CI-інтеграція `check-catalog-fresh`;
  golden-набір складних задач як регрес-тести скілів.

---

## 9. Обслуговування

- Джерело істини незмінне: `docs/mcp/tools/`, `docs/algorithms/`, `MCP_STATUS.md`, `LIMITATIONS.md`.
- Зміна інструмента → перегенерувати каталог → зелений `check-catalog-fresh`.
- Скіли не описують _як влаштований_ алгоритм (це в `docs/algorithms`) — лише _як його застосувати_ до задачі.
