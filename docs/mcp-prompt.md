# MCP Tool Implementation Prompt

## Репозиторій
`d:\Prototype\SuffixTree`

## Чекліст
`docs/mcp-checklist.md` — секція 4 (Per-Tool Checklist). Це єдине джерело правди про те, що реалізовано, а що ні.

## Еталонні файли
- Tool + records: `src/Seqeron/Mcp/Seqeron.Mcp.Parsers/Tools/ParsersTools.cs`
- Тести: `src/Seqeron/Mcp/Seqeron.Mcp.Parsers.Tests/`
- Документація: `docs/mcp/tools/parsers/*.md`, `*.mcp.json`

## Цикл для кожного tool

1. Знайди наступний невиконаний tool у `docs/mcp-checklist.md`
2. Вивчи наявну реалізацію в `Seqeron.Genomics` (DocRef із чекліста)
3. Реалізуй tool + result record за зразком еталонних файлів
4. Створи тести (Schema + Binding)
5. Створи документацію (.md + .mcp.json)
6. `dotnet build && dotnet test` — усі тести мають пройти
7. Познач `[x]` пункти a-j у `docs/mcp-checklist.md`
8. Коміт: `feat(MCP/Parsers): add <tool_name> tool`
9. Повтори з кроку 1

## Критично важливо

- Mission-critical бібліотека. Помилки й неточності неприпустимі.
- Кожен tool має коректно викликати метод із `Seqeron.Genomics`, а не реалізовувати логіку заново.
- Тести мають перевіряти реальну поведінку, а не заглушки.

## Завдання

Почни з наступного невиконаного tool із чекліста.
