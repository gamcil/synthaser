import Select from 'react-select'

/**
 * A domain family filter item in a FilterList component.
 *
 * Props:
 *  @prop data {Object} - object containing data for a rule filter
 *  @prop rule {Object} - the parent rule object
 *  @prop domains {Array} - array of all domains
 *  @prop handleChange {Function} - fn for updating filter in rule
 *  @prop handleRemove {Function} - fn for removing filter from rule
 */
export const FilterItem = props => {
  // Get the domain object corresponding to the selected type, if exists
  const domain = props.domains.find(d => d.name === props.data.type)

  // Generate type options based on domains in the rule
  const typeOptions = props.rule.domains.map(d => ({ label: d, value: d }))

  // Generate domain family options based on families in the selected type
  const domainOptions = domain
    ? domain.domains.map(d => ({ label: d.name, value: d.accession }))
    : []

  const handleChangeName = event => {
    props.handleChange([
      { target: { name: "domains", value: [] } },
      { target: { name: "type", value: event.value } },
    ])
  }
  const handleChangeDomains = event => {
    props.handleChange([
      { target: { name: "domains", value: event ? event.map(d => d.value) : [] } }
    ])
  }

  return (
    <li>
      <button
        type="button"
        onClick={props.handleRemove}
      >
        Delete
      </button>

      <div className="rule-field">
        <label htmlFor="filterName">Domain name:</label>
        <Select
          className="select"
          options={typeOptions}
          onChange={handleChangeName}
          value={domain ? {label: domain.name, value: domain.name} : null}
        />
      </div>

      <div className="rule-field">
        <label htmlFor="filterDomains">Domain types:</label>
        <Select
          className="select"
          options={domainOptions}
          onChange={handleChangeDomains}
          value={props.data.domains.map(d => domainOptions.find(o => o.value === d))}
          isMulti
        />
      </div>
    </li>
  )
}
